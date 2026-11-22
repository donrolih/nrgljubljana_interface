#include <nrgljubljana_interface/solver_core.hpp>

#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <h5/h5.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace nrgljubljana_interface;

TEST(nrgljubljana_interface, siam) { // NOLINT

  // System Parameters
  double U  = 1.0;

  // Construct Parameters
  constr_params_t cp;
  cp.model = "SIAM";
  cp.symtype = "QS";
  cp.mesh_max = 1.0;
  cp.mesh_min = 1e-3;
  cp.mesh_ratio = 1.1;

  // Set up the Solver
  solver_core S(cp);

  // Solve Parameters
  solve_params_t sp;
  sp.T = 1e-3;
  sp.Lambda = 4.0;
  sp.Nz = 2;
  sp.Tmin = 1e-4;
  sp.keep = 50;
  sp.keepenergy = 6.0;

  // Model parameters
  std::map<std::string, double> mp;
  mp["U1"] = 0.5;
  mp["eps1"] = -0.24;
  sp.model_parameters = mp;

  // Low-level NRG Parameters
  nrg_params_t np;
  np.bins = 50;
  S.set_nrg_params(np);

  // Initialize hybridization function
  for (auto w : S.log_mesh)
    S.Delta_w[0][w](0,0) = -0.05i; // a flat band: Gamma = -Im Delta = 0.05

  // Solve the impurity model
  S.solve(sp);

  // Store the Result
  {
    auto arch = h5::file("siam.out.h5", 'w');
    h5_write(arch, "S", S);
  }

  // Compare against the reference data, Gf by Gf, at the precision recorded in siam.ref.json (REF_PRECISION)
  auto out_file = h5::file("siam.out.h5", 'r');
  auto ref_file = h5::file("siam.ref.h5", 'r');
  auto out      = h5::group(out_file).open_group("S");
  auto ref      = h5::group(ref_file).open_group("S");
  for (std::string name : {"A_w", "B_l_w", "B_r_w", "C_w", "Delta_w", "F_l_w", "F_r_w", "G_w", "I_w", "SigmaHartree_w", "Sigma_w",
                           "chi_NN_w", "chi_SS_w"}) {
    ASSERT_EQ(out.has_subgroup(name), ref.has_subgroup(name)) << name << " is in only one of output and reference";
    if (not ref.has_subgroup(name)) continue;
    auto g_out = h5::read<g_w_t>(out, name);
    auto g_ref = h5::read<g_w_t>(ref, name);
    EXPECT_BLOCK_GF_NEAR(g_out, g_ref, REF_PRECISION) << name;
  }
}

MAKE_MAIN
