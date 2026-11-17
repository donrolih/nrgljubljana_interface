// Broadening the stored bins of a result gives what the solve broadened with the same parameters: on the result in
// memory, after an HDF5 round trip, and in a fresh process (DeferredBroadeningFromFile, run by CTest after the
// others, reads only the file they write). S.broaden updates the solver.

#include <nrgljubljana_interface/solver_core.hpp>

#include <triqs/test_tools/gfs.hpp>

#include <h5/h5.hpp>

#include <cmath>
#include <numbers>
#include <string>
#include <vector>

using namespace nrgljubljana_interface;

namespace {

  std::string const result_file = "deferred_broadening.h5";

  std::vector<std::string> const correlators = {"G", "F_l", "F_r", "I", "chi_NN", "chi_SS"};

  solver_core solve_siam() {
    constr_params_t cp;
    cp.model      = "SIAM";
    cp.symtype    = "QS";
    cp.mesh_max   = 1.0;
    cp.mesh_min   = 1e-3;
    cp.mesh_ratio = 1.1;
    solver_core S(cp);

    nrg_params_t np;
    np.bins = 50;
    S.set_nrg_params(np);
    for (auto w : S.log_mesh) S.Delta_w[0][w](0, 0) = 0.05i;

    solve_params_t sp;
    sp.T                = 1e-3;
    sp.Lambda           = 4.0;
    sp.Nz               = 2;
    sp.Tmin             = 1e-4;
    sp.keep             = 50;
    sp.keepenergy       = 6.0;
    sp.model_parameters = {{"U1", 0.5}, {"eps1", -0.24}};
    S.solve(sp);
    if (not S.last_result()) throw std::runtime_error("no result after a solve");
    return S;
  }

  // The solver after a default solve, shared by the tests (S.broaden changes it: copy before)
  solver_core const &solved() {
    static solver_core const S = solve_siam();
    return S;
  }

  double max_diff(g_w_t const &a, g_w_t const &b) {
    double m = 0;
    for (long k = 0; k < a.size(); ++k) m = std::max(m, double(max_element(abs(a[k].data() - b[k].data()))));
    return m;
  }

  // broaden with the parameters of the solve gives the processed Green's functions exactly
  void expect_same_as_solve(result_t const &r) {
    for (auto const &name : correlators) {
      SCOPED_TRACE(name);
      EXPECT_BLOCK_GF_NEAR(broaden(r, name), r.processed.at(name), 0.0);
    }
  }

} // namespace

TEST(DeferredBroadening, SameAsSolve) {
  auto const &r = *solved().last_result();
  expect_same_as_solve(r);
  {
    h5::file f(result_file, 'w');
    h5_write(f, "r", r);
  }
  h5::file f(result_file, 'r');
  auto const s = h5::read<result_t>(f, "r");
  EXPECT_EQ(s.mesh, r.mesh);
  expect_same_as_solve(s);
}

TEST(DeferredBroadening, Parameters) {
  auto const &r = *solved().last_result();
  auto const &G = r.processed.at("G");

  broadening_params_t b;
  b.alpha = r.solve_params.alpha;
  b.gamma = r.solve_params.gamma;
  EXPECT_BLOCK_GF_NEAR(broaden(r, "G", b), G, 0.0);

  // Each parameter reaches the broadening
  b.alpha = 0.6;
  EXPECT_GT(max_diff(broaden(r, "G", b), G), 1e-8);
  b       = {};
  b.gamma = 1.0;
  EXPECT_GT(max_diff(broaden(r, "G", b), G), 1e-8);
  b              = {};
  b.omega0_ratio = 10; // omega0 = 10 T inside the mesh
  EXPECT_GT(max_diff(broaden(r, "G", b), G), 1e-8);

  auto const mesh = triqs::mesh::refreq_log{1e-2, 0.5, 1.2};
  EXPECT_EQ(broaden(r, "G", {}, mesh)[0].mesh(), mesh);
}

TEST(DeferredBroadening, Errors) {
  auto r = *solved().last_result();
  EXPECT_THROW(broaden(r, "X"), std::invalid_argument);
  EXPECT_THROW(broaden(r, "Sigma"), std::invalid_argument);
  EXPECT_THROW(broaden(r, "SigmaHartree"), std::invalid_argument);
  for (auto &t : r.twists) t.spectra.clear();
  EXPECT_THROW(broaden(r, "G"), std::runtime_error);

  constr_params_t cp;
  cp.model   = "SIAM";
  cp.symtype = "QS";
  solver_core S(cp);
  EXPECT_THROW(S.broaden("G"), std::exception); // no result
}

TEST(DeferredBroadening, SolverMethod) {
  solver_core S        = solved();
  auto const G_solve   = *S.G_w;
  auto const F_l_solve = *S.F_l_w;
  auto const Sigma     = *S.Sigma_w;

  broadening_params_t b;
  b.alpha      = 0.6;
  auto const G = S.broaden("G", b);
  EXPECT_BLOCK_GF_NEAR(*S.G_w, G, 0.0);
  EXPECT_BLOCK_GF_NEAR(S.last_result()->processed.at("G"), G, 0.0);
  EXPECT_GT(max_diff(G, G_solve), 1e-8);
  auto A = G;
  for (auto &blk : A.data())
    for (auto &x : blk.data()) x = -x.imag() / std::numbers::pi;
  EXPECT_BLOCK_GF_NEAR(*S.A_w, A, 0.0);
  EXPECT_BLOCK_GF_NEAR(*S.F_l_w, F_l_solve, 0.0); // others keep their broadening
  EXPECT_BLOCK_GF_NEAR(*S.Sigma_w, Sigma, 0.0);

  // The parameters of the solve restore G
  EXPECT_BLOCK_GF_NEAR(S.broaden("G"), G_solve, 0.0);
  EXPECT_BLOCK_GF_NEAR(*S.G_w, G_solve, 0.0);
  EXPECT_THROW(S.broaden("Sigma"), std::exception);
}

// In a fresh process: only the file written by DeferredBroadening.SameAsSolve
TEST(DeferredBroadeningFromFile, SameAsSolve) {
  h5::file f(result_file, 'r');
  expect_same_as_solve(h5::read<result_t>(f, "r"));
}

MAKE_MAIN
