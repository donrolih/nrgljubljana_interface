// Raw results of a solve (solver_core::last_result): what a default solve collects from its twists, that every
// MPI rank holds the same, and what a forward-only and a single-twist solve collect

#include <nrgljubljana_interface/solver_core.hpp>

#include <triqs/test_tools/gfs.hpp>

#include <h5/h5.hpp>

#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <string>

using namespace nrgljubljana_interface;

namespace {

  solver_core solve_siam(options_t const &options = {}) {
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
    S.solve(sp, options);
    if (not S.last_result()) throw std::runtime_error("no result after a solve");
    return S;
  }

  // The result of a default solve, shared by the tests
  result_t const &default_result() {
    static result_t const r = *solve_siam().last_result();
    return r;
  }

  void expect_same_raw(result_t const &a, result_t const &b) {
    EXPECT_EQ(a.z, b.z);
    EXPECT_EQ(a.complete_z_average, b.complete_z_average);
    ASSERT_EQ(a.twists.size(), b.twists.size());
    for (std::size_t k = 0; k < a.twists.size(); ++k) {
      auto const &s = a.twists[k], &t = b.twists[k];
      EXPECT_EQ(s.z, t.z);
      EXPECT_EQ(s.param, t.param);
      EXPECT_EQ(s.spectra, t.spectra);
      EXPECT_EQ(s.custom.columns, t.custom.columns);
      EXPECT_ARRAY_EQ(s.custom.values, t.custom.values);
      EXPECT_EQ(s.td.columns, t.td.columns);
      EXPECT_ARRAY_EQ(s.td.values, t.td.values);
      EXPECT_EQ(s.customfdm, t.customfdm);
      EXPECT_EQ(s.tdfdm, t.tdfdm);
      EXPECT_EQ(s.energy_flow, t.energy_flow);
    }
  }

} // namespace

TEST(RawResults, DefaultSolve) {
  mpi::communicator world;
  auto const &r = default_result();

  EXPECT_EQ(r.z, (std::vector<double>{0.5, 1.0}));
  EXPECT_TRUE(r.complete_z_average);
  EXPECT_EQ(r.model, "SIAM");
  EXPECT_NE(r.manifest.find("[correlators]"), std::string::npos);
  EXPECT_NE(r.nrg_revision.find("nrgljubljana@"), std::string::npos);
  ASSERT_EQ(r.twists.size(), 2u);
  for (auto const &t : r.twists) {
    std::set<std::string> spectra;
    for (auto const &[corr, bins] : t.spectra) {
      spectra.insert(corr);
      EXPECT_FALSE(bins.frequency.empty()) << corr;
      EXPECT_EQ(bins.frequency.size(), bins.weight.size()) << corr;
    }
    EXPECT_EQ(spectra, (std::set<std::string>{"A_d-A_d", "self_d-A_d", "A_d-self_d", "self_d-self_d", "n_d-n_d", "sigma_d-sigma_d"}));
    std::ostringstream z_line; // as in the parameter file
    z_line << std::setprecision(std::numeric_limits<double>::max_digits10) << "\nz=" << t.z << "\n";
    EXPECT_NE(t.param.find(z_line.str()), std::string::npos) << z_line.str();
    EXPECT_GT(t.custom.values.extent(0), 0);
    EXPECT_EQ(t.custom.columns.front(), "T");
    EXPECT_GT(t.td.values.extent(0), 0);
    EXPECT_TRUE(t.customfdm.contains("n_d"));
    EXPECT_TRUE(t.tdfdm.contains("F_fdm"));
    EXPECT_TRUE(t.energy_flow.empty()); // energies = 0
  }
  // The two twists are different discretizations
  EXPECT_NE(r.twists[0].spectra.at("A_d-A_d"), r.twists[1].spectra.at("A_d-A_d"));

  // Every rank holds the same raw results: each writes its own, rank 0 compares
  {
    h5::file f("raw_results_rank" + std::to_string(world.rank()) + ".h5", 'w');
    h5_write(f, "r", r);
  }
  world.barrier();
  if (world.rank() == 0)
    for (int rank = 1; rank < world.size(); ++rank) {
      h5::file f("raw_results_rank" + std::to_string(rank) + ".h5", 'r');
      expect_same_raw(r, h5::read<result_t>(f, "r"));
    }
}

TEST(RawResults, ForwardOnly) {
  options_t options;
  options.fdm      = false;
  options.energies = 5;
  auto const S     = solve_siam(options);
  auto const &r    = *S.last_result();

  ASSERT_EQ(r.twists.size(), 2u);
  for (auto const &t : r.twists) {
    EXPECT_TRUE(t.spectra.empty());
    EXPECT_TRUE(t.customfdm.empty());
    EXPECT_TRUE(t.tdfdm.empty());
    EXPECT_GT(t.custom.values.extent(0), 0);
    ASSERT_FALSE(t.energy_flow.empty());
    for (auto const &step : t.energy_flow) {
      EXPECT_FALSE(step.empty());
      EXPECT_LE(step.size(), 5u);
    }
  }
  EXPECT_TRUE(r.processed.empty());
  EXPECT_TRUE(r.expv.empty());
  EXPECT_FALSE(S.G_w);
  EXPECT_FALSE(S.A_w);
  EXPECT_FALSE(S.Sigma_w);
}

TEST(RawResults, SingleTwist) {
  options_t options;
  options.z     = {0.5};
  auto const S  = solve_siam(options);
  auto const &r = *S.last_result();
  auto const &d = default_result();

  EXPECT_EQ(r.z, (std::vector<double>{0.5}));
  EXPECT_FALSE(r.complete_z_average);
  ASSERT_EQ(r.twists.size(), 1u);
  ASSERT_EQ(d.twists[0].z, 0.5);
  EXPECT_EQ(r.twists[0].spectra, d.twists[0].spectra);
  EXPECT_EQ(r.twists[0].customfdm, d.twists[0].customfdm);
  auto expv = d.twists[0].customfdm; // the average over one twist, without T
  expv.erase("T");
  EXPECT_EQ(r.expv, expv);
  EXPECT_TRUE(r.processed.contains("G"));
  EXPECT_TRUE(r.processed.contains("Sigma"));
}

MAKE_MAIN
