// The impurity solver as the DMFT loop sees it, and the adapter over solver_core.

#pragma once

#include "../common.hpp"
#include "../../solver/options.hpp"
#include "../../solver/result.hpp"
#include "../../templates/manifest.hpp"

#include <mpi/mpi.hpp>

namespace nrgljubljana_interface {
  class solver_core;
}

namespace nrgljubljana_interface::dmft {

  /// The solver as the loop sees it: solve with Δ, and describe the template, the mesh and the ranks
  class impurity_solver_t {
   public:
    virtual ~impurity_solver_t() = default;

    /// Solve with the hybridization Δ; the result is valid until the next solve
    virtual result_t const &solve(g_w_t const &Delta, solve_params_t const &solve_params, options_t const &options) = 0;

    /// The manifest of the template
    [[nodiscard]] virtual templates::manifest_t const &manifest() const = 0;

    /// The mesh of Δ and of the results
    [[nodiscard]] virtual triqs::mesh::refreq_log mesh() const = 0;

    /// The ranks that solve together, and run the loop together
    [[nodiscard]] virtual mpi::communicator communicator() const = 0;
  };

  /// solver_core as an impurity_solver_t; the solver must outlive the adapter
  class solver_core_adapter : public impurity_solver_t {
   public:
    explicit solver_core_adapter(solver_core &solver) : solver_(solver) {}

    result_t const &solve(g_w_t const &Delta, solve_params_t const &solve_params, options_t const &options) override;
    [[nodiscard]] templates::manifest_t const &manifest() const override;
    [[nodiscard]] triqs::mesh::refreq_log mesh() const override;
    [[nodiscard]] mpi::communicator communicator() const override;

   private:
    solver_core &solver_;
  };

} // namespace nrgljubljana_interface::dmft
