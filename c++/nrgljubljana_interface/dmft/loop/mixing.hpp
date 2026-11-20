// Mixing of the hybridization between DMFT iterations. The mixed quantity is the spectral part Γ of Δ: Re Δ follows
// from it and the on-site energy, and is rebuilt after mixing.

#pragma once

#include "../common.hpp"

#include <h5/h5.hpp>

#include <optional>
#include <string>

namespace nrgljubljana_interface::dmft {

  /// What a mixer is given in one iteration
  struct mix_input_t {
    /// The iteration, from 1
    long iteration = 0;
    /// The spectral parts of Δ and of Δ_out from the lattice; same structure and mesh
    g_w_t Gamma, Gamma_out;
    /// μ of the solve, and the filling n(μ), for a mixer that also moves μ
    double mu = 0, filling = 0;
  };

  /// What a mixer did; the same for all mixers
  struct mix_report_t {
    /// "linear", "broyden" or "linear (fallback)"
    std::string method;
    /// ‖Γ_out − Γ‖ = √(Σ_k w_k Σ_ij |ΔΓ_ij(ω_k)|²) over all blocks, w_k the integration weights of the mesh
    double residual_norm = 0;
    /// Broyden: the number of stored iterations, the condition number of its step, whether it fell back to linear
    long history_size = 0;
    double condition  = 0;
    bool fallback     = false;
  };

  /// The result of mixing
  struct mix_result_t {
    /// The mixed Γ, before the causal projection and the rebuild of Δ
    g_w_t Gamma_mix;
    /// μ for the next iteration, if the mixer moves μ
    std::optional<double> mu_next;
    mix_report_t report;
  };

  /// A mixer: the next Γ from the used and the raw one, with a state kept between iterations
  class mixer_t {
   public:
    virtual ~mixer_t() = default;

    virtual mix_result_t mix(mix_input_t const &input) = 0;

    /// Whether the mixer sets μ in this iteration (mix_result_t::mu_next); then the loop takes μ from it
    [[nodiscard]] virtual bool moves_mu(long iteration) const {
      (void)iteration;
      return false;
    }

    /// Write the type, the parameters and the state, for the checkpoint
    virtual void h5_write(h5::group g, std::string const &name) const = 0;

    /// Restore the state; throws std::runtime_error if the group holds another type of mixer
    virtual void h5_read(h5::group g, std::string const &name) = 0;
  };

  /// Linear mixing: Γ_mix = α Γ_out + (1 − α) Γ, element by element over all blocks
  class linear_mixer : public mixer_t {
   public:
    /// Throws std::invalid_argument unless 0 < α ≤ 1
    explicit linear_mixer(double alpha);

    /// Throws std::invalid_argument if Γ and Γ_out differ in structure or mesh, or are not finite
    mix_result_t mix(mix_input_t const &input) override;

    void h5_write(h5::group g, std::string const &name) const override;
    void h5_read(h5::group g, std::string const &name) override;

    [[nodiscard]] double alpha() const { return alpha_; }

   private:
    double alpha_;
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
