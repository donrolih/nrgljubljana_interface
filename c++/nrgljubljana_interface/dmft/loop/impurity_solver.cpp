#include "./impurity_solver.hpp"

#include "../../solver_core.hpp"

#include <stdexcept>

namespace nrgljubljana_interface::dmft {

  result_t const &solver_core_adapter::solve(g_w_t const &Delta, solve_params_t const &solve_params, options_t const &options) {
    solver_.Delta_w = Delta;
    solver_.solve(solve_params, options);
    auto const &r = solver_.last_result();
    if (not r) throw std::runtime_error("solver_core_adapter: the solve left no result");
    return *r;
  }

  templates::manifest_t const &solver_core_adapter::manifest() const { return solver_.template_manifest(); }

  triqs::mesh::refreq_log solver_core_adapter::mesh() const { return solver_.log_mesh; }

  mpi::communicator solver_core_adapter::communicator() const { return solver_.communicator(); }

} // namespace nrgljubljana_interface::dmft
