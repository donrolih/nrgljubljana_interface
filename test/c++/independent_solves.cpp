// Independent solvers on a self communicator: each rank solves its own parameters at the same time, with the same
// results as a serial solve; a solve leaves nrg_params as given, so a second solve repeats the first

#include <nrgljubljana_interface/solver_core.hpp>

#include <triqs/test_tools/gfs.hpp>

#include <h5/h5.hpp>

#include <string>

using namespace nrgljubljana_interface;

namespace {

  double U_of_rank(int rank) { return 0.5 + 0.1 * rank; }

  solver_core make_solver(mpi::communicator c) {
    constr_params_t cp;
    cp.model      = "SIAM";
    cp.symtype    = "QS";
    cp.mesh_max   = 1.0;
    cp.mesh_min   = 1e-3;
    cp.mesh_ratio = 1.1;
    solver_core S(cp, c);
    nrg_params_t np;
    np.bins = 50;
    S.set_nrg_params(np);
    for (auto w : S.log_mesh) S.Delta_w[0][w](0, 0) = 0.05i;
    return S;
  }

  solve_params_t solve_params(double U) {
    solve_params_t sp;
    sp.T                = 1e-3;
    sp.Lambda           = 4.0;
    sp.Nz               = 2;
    sp.Tmin             = 1e-4;
    sp.keep             = 50;
    sp.keepenergy       = 6.0;
    sp.model_parameters = {{"U1", U}, {"eps1", -U / 2}};
    return sp;
  }

  result_t solve(solver_core &S, double U) {
    S.solve(solve_params(U));
    return *S.last_result();
  }

  void expect_same_raw(result_t const &a, result_t const &b) {
    EXPECT_EQ(a.z, b.z);
    ASSERT_EQ(a.twists.size(), b.twists.size());
    for (std::size_t k = 0; k < a.twists.size(); ++k) {
      auto const &s = a.twists[k], &t = b.twists[k];
      EXPECT_EQ(s.param, t.param);
      EXPECT_EQ(s.spectra, t.spectra);
      EXPECT_ARRAY_EQ(s.custom.values, t.custom.values);
      EXPECT_ARRAY_EQ(s.td.values, t.td.values);
      EXPECT_EQ(s.customfdm, t.customfdm);
      EXPECT_EQ(s.tdfdm, t.tdfdm);
    }
  }

} // namespace

TEST(IndependentSolves, SelfCommunicator) {
  mpi::communicator world;
  auto S = make_solver(mpi::communicator{MPI_COMM_SELF});
  auto const U = U_of_rank(world.rank());

  // Every rank solves its own U, at the same time and in the same directory
  auto const np_given = S.nrg_params;
  auto const first    = solve(S, U);
  EXPECT_EQ(S.nrg_params.Nmax, np_given.Nmax);
  EXPECT_EQ(S.nrg_params.mMAX, np_given.mMAX);
  EXPECT_EQ(S.nrg_params.xmax, np_given.xmax);
  EXPECT_EQ(S.last_solve_params->bandrescale, 1.0); // resolved from mesh_max

  // A second solve on the same solver repeats the first
  expect_same_raw(first, solve(S, U));

  // A serial solve of another rank's U gives that rank's results
  {
    h5::file f("independent_solves_rank" + std::to_string(world.rank()) + ".h5", 'w');
    h5_write(f, "r", first);
  }
  world.barrier();
  if (world.rank() == 0)
    for (int rank = 1; rank < world.size(); ++rank) {
      auto serial = make_solver(mpi::communicator{MPI_COMM_SELF});
      h5::file f("independent_solves_rank" + std::to_string(rank) + ".h5", 'r');
      expect_same_raw(solve(serial, U_of_rank(rank)), h5::read<result_t>(f, "r"));
    }
}

MAKE_MAIN
