// The chemical potential between DMFT iterations: fixed, or adjusted towards a filling with the self-energy held
// fixed, as dmft-jamova's occupancy_mode=accurate.

#pragma once

#include "../common.hpp"
#include "../lattice/lattice.hpp"

#include <h5/h5.hpp>

#include <map>
#include <optional>
#include <string>

namespace nrgljubljana_interface::dmft {

  /// Parameters of mu_control_t
  struct mu_control_params_t {
    /// The target filling (the total of occupancy); unset: μ is kept fixed
    std::optional<double> goal = {};
    /// The largest step of μ in one iteration
    double maxdx = 0.05;
    /// The damping of the step to the root, 0 < under ≤ 1; the bracket is searched up to maxdx / under
    double under = 1.0;
    /// The allowed |n − goal| of the root, and the accuracy of the root in μ
    double root_tol = 1e-6, mu_tol = 1e-7;
    /// Lattice evaluations per iteration, at least 4
    long max_evaluations = 12;
    /// The allowed |∫ A_ii dω − 1| of the lattice spectral function
    double weight_tol = 1e-4;
    /// Tolerances of the Fermi integral of the filling
    double epsabs = 1e-10, epsrel = 1e-9;
  };

  void h5_write(h5::group g, std::string const &name, mu_control_params_t const &p);
  void h5_read(h5::group g, std::string const &name, mu_control_params_t &p);

  /// What the μ control is given in one iteration
  struct mu_input_t {
    /// The self-energy, held fixed
    g_w_t Sigma;
    /// μ of the solve
    double mu_used = 0;
    /// The temperature, and the degeneracy of each block (the manifest's)
    double T = 0;
    std::map<std::string, long> degeneracy;
  };

  /// What the μ control did
  struct mu_report_t {
    /// "fixed", "within_tolerance", "bracketed" or "limited"
    std::string status;
    /// The filling at μ_used and at μ_next
    double filling_used = 0, filling_next = 0;
    /// The largest |∫ A_ii dω − 1| at μ_used and at μ_next
    double weight_error_used = 0, weight_error_next = 0;
    /// The root of n(μ) = goal, when bracketed
    std::optional<double> mu_root;
    /// The number of lattice evaluations
    long evaluations = 0;
  };

  /// The result of adjusting μ
  struct mu_result_t {
    double mu_next = 0;
    /// The filling at μ_used
    double filling = 0;
    mu_report_t report;
  };

  /**
   * The chemical potential: fixed without a goal, otherwise adjusted towards the goal with Σ held fixed.
   *
   * With a goal: the filling at μ_used; within root_tol of the goal μ is kept. Otherwise the filling at
   * μ_used ± maxdx / under in the direction of the error, which must not move the filling the wrong way. If the two
   * bracket the goal, Brent's method (GSL) to mu_tol, the residual checked against root_tol, and
   * μ_next = μ_used + clamp(under (μ* − μ_used), ±maxdx); otherwise μ_next = μ_used ± maxdx. At μ_next the filling
   * error must not have grown by more than root_tol. Evaluations are cached and limited to max_evaluations.
   *
   * With a goal every evaluation checks that each ∫ A_ii dω is 1 within weight_tol and that the filling lies between
   * 0 and the number of states; without one the weight error is only reported. Errors: std::invalid_argument for
   * invalid parameters or input, std::runtime_error for a failed search.
   */
  class mu_control_t {
   public:
    explicit mu_control_t(mu_control_params_t const &params = {});

    [[nodiscard]] mu_result_t adjust(mu_input_t const &input, lattice_t const &lattice) const;

    [[nodiscard]] mu_control_params_t const &params() const { return params_; }

   private:
    mu_control_params_t params_;
  };

} // namespace nrgljubljana_interface::dmft
