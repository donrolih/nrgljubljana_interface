// The DMFT loop: the impurity solver, the lattice and the loop controls (mixer, μ control, convergence test).

#pragma once

#include "../causality.hpp"
#include "../common.hpp"
#include "../lattice/lattice.hpp"
#include "./convergence.hpp"
#include "./impurity_solver.hpp"
#include "./mixing.hpp"
#include "./mu_control.hpp"

#include <triqs/utility/macros.hpp>

#include <complex>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <string>

namespace nrgljubljana_interface::dmft {

  /// Parameters of the loop
  struct loop_params_t {
    /// Solve parameters (T, Λ, Nz, the model parameters except the impurity levels)
    solve_params_t solve_params = {};
    /// What each solve computes; must give G, F_l, F_r and I broadened, for the full set of z
    options_t options = {};
    /// Floors of the spectral parts of Δ and Σ
    double Gamma_floor = 1e-6, Sigma_floor = 1e-12;
    /// Directory of the files of the run (rank 0): checkpoint.h5, history.h5, stats.dat; empty: no files
    std::string run_dir = {};
    /// Append each iteration to history.h5
    bool keep_history = true;

    static std::string hdf5_format() { return "NRGLJUBLJANA_INTERFACE_LoopParams"; }
  };

  void h5_write(h5::group g, std::string const &name, loop_params_t const &p);
  void h5_read(h5::group g, std::string const &name, loop_params_t &p);

  /// One iteration of the loop
  struct step_result_t {
    long iteration = 0;
    /// μ of the solve, μ for the next iteration, and the filling at μ
    double mu = 0, mu_next = 0, filling = 0;
    /// The causal projections of Σ and of the next Δ
    causality_report_t Sigma_projection = {}, Delta_projection = {};
    mu_report_t mu_report = {};
    mix_report_t mix_report = {};
    convergence_result_t convergence = {};
    /// Expectation values of the solve
    std::map<std::string, double> expv = {};
    /// The next Δ could not be prepared after a terminal status; Δ is kept
    bool kept_Delta = false;

    static std::string hdf5_format() { return "NRGLJUBLJANA_INTERFACE_StepResult"; }
  };

  void h5_write(h5::group g, std::string const &name, step_result_t const &r);
  void h5_read(h5::group g, std::string const &name, step_result_t &r);

  /**
   * The DMFT loop. One step: the impurity levels ⟨h⟩_ii − μ from the manifest's [levels]; the solve with Δ; Σ
   * made causal; G_loc at μ; μ_next from the μ control (or from the mixer, when it moves μ); (G_loc, Δ_out) at
   * μ_next; the convergence status with G_imp = impurity_gf(Δ, Σ, −μ); Γ and Γ_out mixed and the next Δ
   * from hybridization_from_Gamma. Rank 0 broadcasts Δ_mix, μ_next and the status.
   *
   * Without a μ control μ stays fixed and the filling is still measured. Construction throws std::invalid_argument
   * for options that do not give G, F_l, F_r and I broadened for the full set of z, a manifest without a level for
   * every diagonal element or with Δ and G of different structure, a level also given as a model parameter, and a
   * lattice that cannot give the on-site energy of the structure.
   */
  class dmft_loop {
   public:
    using callback_t = std::function<void(step_result_t const &)>;

    C2PY_IGNORE dmft_loop(impurity_solver_t &solver, lattice_t const &lattice, mixer_t &mixer, convergence_t &convergence, loop_params_t params,
                          std::optional<mu_control_t> const &mu_control = {});

    /// With the solver of the package, through a solver_core_adapter held by the loop
    dmft_loop(solver_core &solver, lattice_t const &lattice, mixer_t &mixer, convergence_t &convergence, loop_params_t params,
              std::optional<mu_control_t> const &mu_control = {});

    /// Δ from the lattice at Σ0 and μ0, through hybridization_from_Gamma; Σ0 in the structure and mesh of G. The loop
    /// controls are reset, and a history.h5 in the run directory is removed.
    void initialize(g_w_t const &Sigma0, double mu0);

    /// As initialize, with Σ0 times the unit matrix in every block and at every frequency
    void initialize(std::complex<double> Sigma0, double mu0);

    /**
     * Continue the run whose checkpoint.h5 is in run_dir: the next step is the iteration after it, with its Δ_mix and
     * μ_next, and the states of the mixer and the convergence test restored. The stored parameters must equal the
     * current ones (solve parameters, options, floors, μ-control and convergence parameters, the descriptions of the
     * mixer and the lattice, the manifest, the mesh); otherwise std::invalid_argument names each difference, unless
     * allow_changes, which continues with the current parameters.
     */
    void resume(std::string const &run_dir, bool allow_changes = false);

    /// Start a new run from the Δ_mix and μ_next of the checkpoint in run_dir: iteration 0, the loop controls reset
    void initialize_from(std::string const &run_dir);

    /// One iteration; throws std::logic_error before initialize
    step_result_t step();

    /// Steps until the status is converged or stop, or max_steps (0: no limit) steps were taken; the last status
    std::string run(long max_steps = 0);

    /// Called with the result of each step, after it is committed
    C2PY_IGNORE void set_callback(callback_t callback) { callback_ = std::move(callback); }

    [[nodiscard]] long iteration() const { return iteration_; }
    [[nodiscard]] double mu() const { return mu_; }
    [[nodiscard]] C2PY_IGNORE g_w_t const &Delta() const;
    [[nodiscard]] C2PY_IGNORE std::optional<g_w_t> const &Sigma() const { return Sigma_; }
    [[nodiscard]] C2PY_IGNORE std::optional<g_w_t> const &G_loc() const { return G_loc_; }
    [[nodiscard]] C2PY_IGNORE std::optional<step_result_t> const &last_step() const { return last_; }
    [[nodiscard]] block_matrices_t const &onsite_energy() const { return onsite_; }
    [[nodiscard]] loop_params_t const &params() const { return params_; }

   private:
    std::unique_ptr<impurity_solver_t> owned_solver_;
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
    std::optional<g_w_t> Delta_, Sigma_, G_loc_, G_loc_previous_;
    std::optional<step_result_t> last_;

    void setup();
    [[nodiscard]] solve_params_t with_levels(double mu) const;
    void write_statistics(step_result_t const &r) const;
    void write_checkpoint(step_result_t const &r, result_t const &result, g_w_t const &Delta_mix, g_w_t const &Sigma, g_w_t const &G_loc,
                          g_w_t const &G_imp) const;
    void append_history(step_result_t const &r, g_w_t const &Delta_out, g_w_t const &Sigma, g_w_t const &G_loc) const;
    void start_files() const;
  };

} // namespace nrgljubljana_interface::dmft
