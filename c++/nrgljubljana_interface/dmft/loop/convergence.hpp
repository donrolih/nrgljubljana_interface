// Convergence of the DMFT loop: continue, converged or stop, as dmft-jamova's convergence test.

#pragma once

#include "../common.hpp"

#include <h5/h5.hpp>
#include <triqs/utility/macros.hpp>

#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace nrgljubljana_interface::dmft {

  /// Parameters of convergence_t
  struct convergence_params_t {
    /// Every distance of G_loc between consecutive iterations in the window must be below eps
    double eps = 1e-6;
    /// The number of consecutive distances in the window
    long window = 3;
    /// The first and the last iteration at which the loop may converge; at max_iter it stops
    long min_iter = 1, max_iter = 50;
    /// With a goal: the allowed |n(μ) − goal| and |μ_next − μ|
    double filling_tol = 1e-6, mu_tol = 1e-7;
    /// When set, the distance between G_imp and G_loc must also be at most this
    std::optional<double> eps_imp_loc = {};

    bool operator==(convergence_params_t const &) const = default;
    static std::string hdf5_format() { return "NRGLJUBLJANA_INTERFACE_ConvergenceParams"; }
  };

  void h5_write(h5::group g, std::string const &name, convergence_params_t const &p);
  void h5_read(h5::group g, std::string const &name, convergence_params_t &p);

  /// What the convergence test is given in one iteration
  struct convergence_input_t {
    /// The iteration, from 1; consecutive
    long iteration = 0;
    /// G_loc at μ, and that of the previous iteration (unset in the first)
    g_w_t G_loc;
    std::optional<g_w_t> G_loc_previous = {};
    /// The impurity Green's function
    g_w_t G_imp;
    /// The filling at μ; the goal when μ is controlled; μ_next − μ
    double filling = 0;
    std::optional<double> goal = {};
    double mu_step = 0;
  };

  /// The status and the measures behind it
  struct convergence_result_t {
    /// "continue", "converged" or "stop"
    std::string status = {};
    /// The largest spectral distance of G_loc to the previous iteration (unset in the first)
    std::optional<double> distance = {};
    /// The largest distance in the window, and the number of distances in it
    double window_max   = 0;
    long window_samples = 0;
    /// The largest spectral distance between G_imp and G_loc
    double distance_imp_loc = 0;
    /// Which criteria hold (an inactive criterion holds)
    bool window_ok = false, iterations_ok = false, filling_ok = false, mu_ok = false, imp_loc_ok = false;

    static std::string hdf5_format() { return "NRGLJUBLJANA_INTERFACE_ConvergenceResult"; }
  };

  void h5_write(h5::group g, std::string const &name, convergence_result_t const &r);
  void h5_read(h5::group g, std::string const &name, convergence_result_t &r);

  /**
   * The convergence test. converged: a full window of distances all strictly below eps, min_iter ≤ n ≤ max_iter,
   * with a goal |n − goal| ≤ filling_tol and |μ_next − μ| ≤ mu_tol, and with eps_imp_loc the G_imp–G_loc
   * distance at most that; stop: not converged at max_iter; continue otherwise. The state is the history of distances.
   *
   * Throws std::invalid_argument for invalid parameters, a non-consecutive iteration, a missing previous G_loc after
   * the first iteration, and Green's functions that differ in structure or mesh.
   */
  class convergence_t {
   public:
    explicit convergence_t(convergence_params_t const &params = {});

    convergence_result_t decide(convergence_input_t const &input);

    /// Forget the history, for a new run
    void reset() {
      history_.clear();
      last_iteration_ = 0;
    }

    /// The parameters and the history of distances, for the checkpoint
    C2PY_IGNORE void h5_write(h5::group g, std::string const &name) const;
    C2PY_IGNORE void h5_read(h5::group g, std::string const &name);

    [[nodiscard]] convergence_params_t const &params() const { return params_; }
    /// (iteration, distance) of every iteration after the first
    [[nodiscard]] std::vector<std::pair<long, double>> const &history() const { return history_; }

   private:
    convergence_params_t params_;
    std::vector<std::pair<long, double>> history_;
    long last_iteration_ = 0;
  };

} // namespace nrgljubljana_interface::dmft
