// The DMFT loop: the impurity solver, the lattice and the loop controls (mixer, μ control, convergence test).

#pragma once

#include "../causality.hpp"
#include "../common.hpp"
#include "../lattice/lattice.hpp"
#include "./convergence.hpp"
#include "./impurity_solver.hpp"
#include "./mixing.hpp"
#include "./mu_control.hpp"

#include <functional>
#include <map>
#include <optional>
#include <string>

namespace nrgljubljana_interface::dmft {

  /// Parameters of the loop
  struct loop_params_t {
    /// Solve parameters (T, Λ, Nz, the model parameters except the impurity levels)
    solve_params_t solve_params;
    /// What each solve computes; must give G, F_l, F_r and I broadened, for the full set of z
    options_t options;
    /// Floors of the spectral parts of Δ and Σ
    double gamma_floor = 1e-6, sigma_floor = 1e-12;
    /// Directory of the statistics table (rank 0); empty: no files
    std::string run_dir;
  };

  /// One iteration of the loop
  struct step_result_t {
    long iteration = 0;
    /// μ of the solve, μ for the next iteration, and the filling at μ_used
    double mu_used = 0, mu_next = 0, filling = 0;
    /// The causal projections of Σ and of the next Δ
    causality_report_t sigma_projection, delta_projection;
    mu_report_t mu;
    mix_report_t mix;
    convergence_result_t convergence;
    /// Expectation values of the solve
    std::map<std::string, double> expv;
    /// The next Δ could not be prepared after a terminal decision; Δ_used is kept
    bool kept_delta = false;
  };

  /**
   * The DMFT loop. One step: the impurity levels ⟨h⟩_ii − μ from the manifest's [levels]; the solve with Δ_used; Σ
   * made causal; G_loc at μ_used; μ_next from the μ control (or from the mixer, when it moves μ); (G_loc, Δ_raw) at
   * μ_next; the convergence decision with G_imp = impurity_gf(Δ_used, Σ, −μ); Γ_used and Γ_raw mixed and the next Δ
   * from hybridization_from_gamma. Rank 0 broadcasts Δ_next, μ_next and the decision.
   *
   * Without a μ control μ stays fixed and the filling is still measured. Construction throws std::invalid_argument
   * for options that do not give G, F_l, F_r and I broadened for the full set of z, a manifest without a level for
   * every diagonal element or with Δ and G of different structure, a level also given as a model parameter, and a
   * lattice that cannot give the on-site energy of the structure.
   */
  class dmft_loop {
   public:
    using callback_t = std::function<void(dmft_loop const &, step_result_t const &)>;

    dmft_loop(impurity_solver_t &solver, lattice_t const &lattice, mixer_t &mixer, convergence_t &convergence, loop_params_t params,
              std::optional<mu_control_t> const &mu_control = {});

    /// Δ_used from the lattice at Σ0 and μ0, through hybridization_from_gamma; Σ0 in the structure and mesh of G
    void initialize(g_w_t const &Sigma0, double mu0);

    /// One iteration; throws std::logic_error before initialize
    step_result_t step();

    /// Steps until the decision is converged or stop, or max_steps (0: no limit) steps were taken; the last decision
    std::string run(long max_steps = 0);

    void set_callback(callback_t callback) { callback_ = std::move(callback); }

    [[nodiscard]] long iteration() const { return iteration_; }
    [[nodiscard]] double mu() const { return mu_; }
    [[nodiscard]] g_w_t const &Delta_used() const;
    [[nodiscard]] std::optional<g_w_t> const &Sigma() const { return Sigma_; }
    [[nodiscard]] std::optional<g_w_t> const &G_loc() const { return G_loc_; }
    [[nodiscard]] std::optional<step_result_t> const &last_step() const { return last_; }
    [[nodiscard]] block_matrices_t const &onsite_energy() const { return onsite_; }
    [[nodiscard]] loop_params_t const &params() const { return params_; }

   private:
    impurity_solver_t &solver_;
    lattice_t const &lattice_;
    mixer_t &mixer_;
    convergence_t &convergence_;
    loop_params_t params_;
    mu_control_t mu_control_;
    block_matrices_t onsite_;
    std::map<std::string, long> degeneracy_;
    callback_t callback_;

    long iteration_ = 0;
    double mu_ = 0, last_mu_step_ = 0;
    std::optional<g_w_t> Delta_used_, Sigma_, G_loc_, G_loc_previous_;
    std::optional<step_result_t> last_;

    [[nodiscard]] solve_params_t with_levels(double mu) const;
    void write_statistics(step_result_t const &r) const;
  };

} // namespace nrgljubljana_interface::dmft
