// Mixing of the hybridization between DMFT iterations. The mixed quantity is the spectral part Γ of Δ: Re Δ follows
// from it and the on-site energy, and is rebuilt after mixing.

#pragma once

#include "../causality.hpp"
#include "../common.hpp"

#include <h5/h5.hpp>
#include <triqs/utility/macros.hpp>

#include <optional>
#include <string>
#include <vector>

namespace nrgljubljana_interface::dmft {

  /// What a mixer is given in one iteration
  struct mix_input_t {
    /// The iteration, from 1
    long iteration = 0;
    /// The spectral parts of Δ and of Δ_out from the lattice; same structure and mesh
    g_w_t Gamma, Gamma_out;
    /// μ of the solve, and the filling n(μ), for a mixer that also moves μ
    double mu = 0, filling = 0;
    /// The target filling, for a mixer that also moves μ
    std::optional<double> goal = {};
  };

  /// What a mixer did; the same for all mixers
  struct mix_report_t {
    /// "linear", "broyden" or "linear (fallback)"
    std::string method = {};
    /// ‖Γ_out − Γ‖ = √(Σ_k w_k Σ_ij |ΔΓ_ij(ω_k)|²) over all blocks, w_k the integration weights of the mesh
    double residual_norm = 0;
    /// Broyden: the number of stored iterations, the condition number of its step, whether it fell back to linear
    long history_size = 0;
    double condition  = 0;
    bool fallback     = false;

    static std::string hdf5_format() { return "NRGLJUBLJANA_INTERFACE_MixReport"; }
  };

  void h5_write(h5::group g, std::string const &name, mix_report_t const &r);
  void h5_read(h5::group g, std::string const &name, mix_report_t &r);

  /// The result of mixing
  struct mix_result_t {
    /// The mixed Γ, before the causal projection and the rebuild of Δ
    g_w_t Gamma_mix;
    /// μ for the next iteration, if the mixer moves μ
    std::optional<double> mu_next = {};
    mix_report_t report = {};
  };

  /// The result of mixing hybridizations
  struct hybridization_mix_t {
    /// The next Δ: the mixed Γ with the floor, and the Hermitian part rebuilt from it and the on-site energy
    g_w_t Delta_mix;
    /// μ for the next iteration, if the mixer moves μ
    std::optional<double> mu_next = {};
    mix_report_t report = {};
    /// The causal projection of the mixed Γ
    causality_report_t projection = {};
  };

  /// A mixer: the next Γ from the used and the raw one, with a state kept between iterations
  class mixer_t {
   public:
    virtual ~mixer_t() = default;

    virtual mix_result_t mix(mix_input_t const &input) = 0;

    /**
     * Mix two hybridizations: Γ and Γ_out, the spectral parts of Δ and Δ_out, mixed by mix, and the next Δ from
     * hybridization_from_Gamma(Γ_mix, onsite_energy, Gamma_floor). μ, the filling and the goal are those of
     * mix_input_t, for a mixer that moves μ. Errors as for mix and hybridization_from_Gamma.
     */
    hybridization_mix_t mix_hybridization(long iteration, g_w_t const &Delta, g_w_t const &Delta_out, block_matrices_t const &onsite_energy,
                                          double Gamma_floor, double mu = 0, double filling = 0, std::optional<double> goal = {});

    /// Whether the mixer sets μ in this iteration (mix_result_t::mu_next); then the loop takes μ from it
    [[nodiscard]] virtual bool moves_mu(long iteration) const {
      (void)iteration;
      return false;
    }

    /// Forget the history, for a new run
    virtual void reset() {}

    /// The type and the parameters, e.g. to check on resume that the mixer is the same
    [[nodiscard]] virtual std::string description() const = 0;

    /// Write the type, the parameters and the state, for the checkpoint
    C2PY_IGNORE virtual void h5_write(h5::group g, std::string const &name) const = 0;

    /// Restore the state; throws std::runtime_error if the group holds another type of mixer
    C2PY_IGNORE virtual void h5_read(h5::group g, std::string const &name) = 0;
  };

  /// Linear mixing: Γ_mix = α Γ_out + (1 − α) Γ, element by element over all blocks
  class linear_mixer : public mixer_t {
   public:
    /// Throws std::invalid_argument unless 0 < α ≤ 1
    explicit linear_mixer(double alpha);

    /// Throws std::invalid_argument if Γ and Γ_out differ in structure or mesh, or are not finite
    mix_result_t mix(mix_input_t const &input) override;

    [[nodiscard]] std::string description() const override;
    C2PY_IGNORE void h5_write(h5::group g, std::string const &name) const override;
    C2PY_IGNORE void h5_read(h5::group g, std::string const &name) override;

    [[nodiscard]] double alpha() const { return alpha_; }

   private:
    double alpha_;
  };

  /// Parameters of broyden_mixer
  struct broyden_params_t {
    /// The initial inverse Jacobian of Γ (the linear mixing factor), 0 < alpha ≤ 1
    double alpha = 0.5;
    /// The number of iterations kept, at least 1
    long history = 10;
    /// Include μ in the vector, after warmup iterations; needs a μ control with a goal
    bool mu = false;
    /// The initial inverse Jacobian of μ (the step per unit of goal − n), positive
    double mu_alpha = 0.5;
    /// The iterations in which the μ control moves μ, before the mixer does
    long warmup = 0;
    /// The largest step of μ in one iteration
    double max_mu_step = 0.05;
    /// Weight the components of Γ by √w_k, w_k the integration weights of the mesh; false: unweighted
    bool weighted = true;
    /// The regularization added to the diagonal of the normalized least-squares system
    double w0 = 1e-4;
    /// A step whose system has a larger condition number falls back to linear mixing
    double max_condition = 1e12;

    bool operator==(broyden_params_t const &) const = default;
    static std::string hdf5_format() { return "NRGLJUBLJANA_INTERFACE_BroydenParams"; }
  };

  void h5_write(h5::group g, std::string const &name, broyden_params_t const &p);
  void h5_read(h5::group g, std::string const &name, broyden_params_t &p);

  /**
   * Broyden mixing (the modified, multisecant form). The vector x holds the independent real components of Γ — the
   * diagonal, the real and imaginary parts of the upper triangle — of all blocks and mesh points, each times √w_k
   * when weighted, and with mu after warmup also μ; the residual F is Γ_out − Γ, and goal − n for μ. Over the last
   * history iterations, with ΔF_i the normalized differences of consecutive residuals and the factors a (alpha for
   * Γ, mu_alpha for μ):
   *
   *   x_next = x + a F − Σ_i (a ΔF_i + Δx_i) c_i,   (ΔFᵀΔF + w0) c = ΔFᵀ F.
   *
   * Pairs with a vanishing ΔF are skipped. The step falls back to the linear one, x + a F, if it is not finite, if
   * the condition number of the system exceeds max_condition, or if an eigenvalue of Γ_mix at a mesh point is below
   * −1e−12 max(1, max |Γ|). The μ step is clamped to max_mu_step. The history holds the iterations from the first,
   * with μ and n, and must be consecutive. Errors: std::invalid_argument for invalid parameters or input.
   */
  class broyden_mixer : public mixer_t {
   public:
    explicit broyden_mixer(broyden_params_t const &params = {});

    mix_result_t mix(mix_input_t const &input) override;
    [[nodiscard]] bool moves_mu(long iteration) const override { return params_.mu and iteration > params_.warmup; }
    void reset() override { history_.clear(); }
    [[nodiscard]] std::string description() const override;
    C2PY_IGNORE void h5_write(h5::group g, std::string const &name) const override;
    C2PY_IGNORE void h5_read(h5::group g, std::string const &name) override;

    [[nodiscard]] broyden_params_t const &params() const { return params_; }
    /// The number of iterations kept
    [[nodiscard]] long history_size() const { return long(history_.size()); }

   private:
    struct sample_t {
      long iteration = 0;
      std::vector<double> x, F; // the components of Γ, weighted
      double mu = 0, filling = 0;
    };
    broyden_params_t params_;
    std::vector<sample_t> history_;
  };

  namespace detail {

    /// The integration weights of the mesh points: half the distance to each neighbour
    std::vector<double> integration_weights(std::vector<double> const &points);

    /// √(Σ_k w_k Σ_ij |a_ij(ω_k) − b_ij(ω_k)|²) over all blocks
    double weighted_distance(g_w_t const &a, g_w_t const &b);

    /// Throws std::runtime_error unless the group was written by a mixer of the given type
    void check_mixer_type(h5::group g, std::string const &name, std::string const &type);

  } // namespace detail

} // namespace nrgljubljana_interface::dmft
