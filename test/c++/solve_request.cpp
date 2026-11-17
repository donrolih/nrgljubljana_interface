// Solve requests and results (solver/request.hpp, solver/result.hpp): validation and HDF5 round trips

#include <nrgljubljana_interface/solver/request.hpp>
#include <nrgljubljana_interface/solver/result.hpp>

#include <triqs/test_tools/gfs.hpp>

#include <h5/h5.hpp>

#include <stdexcept>
#include <string>
#include <vector>

using namespace nrgljubljana_interface;

namespace {

  templates::manifest_t const manifest = templates::parse_manifest(R"([template]
format=1
model=SIAM
symtype=QS
channels=1
ops=n_d A_d self_d sigma_d SigmaHartree
params=U1 eps1
[structure]
gf=imp:1
degeneracy=imp:2
[bath]
coefficient_channels=imp[0,0]
[correlators]
G imp[0,0]=A_d-A_d
F_l imp[0,0]=self_d-A_d
F_r imp[0,0]=A_d-self_d
I imp[0,0]=self_d-self_d
chi_NN imp[0,0]=n_d-n_d
[hartree]
imp[0,0]=SigmaHartree
)");

  // The error message of validate_request, or "" if the request is valid
  std::string error_of(solve_request_t const &r) {
    try {
      validate_request(r, manifest);
    } catch (std::invalid_argument const &e) { return e.what(); }
    return {};
  }

  bool contains(std::string const &s, std::string const &part) { return s.find(part) != std::string::npos; }

  using names = std::vector<std::string>;

} // namespace

// ---- Validation ----

TEST(SolveRequest, DefaultsAreValid) {
  EXPECT_EQ(error_of({}), "");
  EXPECT_EQ(requested_correlators({}, manifest), (names{"F_l", "F_r", "G", "I", "chi_NN"}));
}

TEST(SolveRequest, UnsetCorrelatorsFollowFdm) {
  solve_request_t r;
  r.fdm = false;
  EXPECT_EQ(error_of(r), "");
  EXPECT_TRUE(requested_correlators(r, manifest).empty());
  r.correlators = names{};
  EXPECT_EQ(error_of(r), "");
}

TEST(SolveRequest, CorrelatorsNeedFdm) {
  solve_request_t r;
  r.fdm         = false;
  r.correlators = names{"G"};
  EXPECT_TRUE(contains(error_of(r), "correlators require fdm=true"));
}

TEST(SolveRequest, Correlators) {
  solve_request_t r;
  r.correlators = names{"G", "chi_NN"};
  EXPECT_EQ(error_of(r), "");
  EXPECT_EQ(requested_correlators(r, manifest), (names{"G", "chi_NN"}));
  r.correlators = names{"G", "chi_SS"};
  EXPECT_TRUE(contains(error_of(r), "unknown correlator \"chi_SS\""));
  r.correlators = names{"G", "G"};
  EXPECT_TRUE(contains(error_of(r), "correlator \"G\" given twice"));
  r.correlators = names{"A_d-A_d"}; // NRG correlator strings are not logical names
  EXPECT_TRUE(contains(error_of(r), "unknown correlator"));
}

TEST(SolveRequest, BroadenAndBins) {
  solve_request_t r;
  r.broaden = false;
  EXPECT_EQ(error_of(r), "");
  r.keep_bins = false;
  EXPECT_TRUE(contains(error_of(r), "broaden=false requires keep_bins=true"));
  r.broaden = true;
  EXPECT_EQ(error_of(r), "");
}

TEST(SolveRequest, EnergiesAndZ) {
  solve_request_t r;
  r.energies = -1;
  EXPECT_TRUE(contains(error_of(r), "energies must not be negative"));
  r.energies = 10;
  r.z        = {1.0};
  EXPECT_EQ(error_of(r), "");
  r.z = {0.25, 0.5, 1.0};
  EXPECT_EQ(error_of(r), "");
  for (double z : {0.0, -0.5, 1.5}) {
    r.z = {0.5, z};
    EXPECT_TRUE(contains(error_of(r), "is outside (0, 1]")) << z;
  }
  r.z = {0.5, 0.5};
  EXPECT_TRUE(contains(error_of(r), "given twice"));
}

// ---- HDF5 ----

TEST(SolveRequest, H5RoundTrip) {
  solve_request_t r;
  r.fdm          = true;
  r.correlators  = names{"G", "chi_NN"};
  r.broaden      = false;
  r.energies     = 7;
  r.z            = {0.5, 1.0};
  r.keep_workdir = true;
  {
    h5::file f("solve_request.h5", 'w');
    h5_write(f, "r", r);
    h5_write(f, "default", solve_request_t{});
    h5_write(f, "none", solve_request_t{.correlators = names{}});
  }
  h5::file f("solve_request.h5", 'r');
  EXPECT_EQ(h5::read<solve_request_t>(f, "r"), r);
  EXPECT_EQ(h5::read<solve_request_t>(f, "default"), solve_request_t{}); // unset correlators stay unset
  EXPECT_EQ(h5::read<solve_request_t>(f, "none"), solve_request_t{.correlators = names{}});
}

namespace {

  solve_result_t sample_result() {
    solve_result_t r;
    r.request.correlators = names{"G"};
    r.request.energies    = 2;
    r.solve_params.T      = 1e-3;
    r.solve_params.Nz     = 2;
    r.solve_params.model_parameters = {{"U1", 0.5}, {"eps1", -0.25}};
    r.nrg_params.bins     = 50;
    r.model               = "SIAM";
    r.symtype             = "QS";
    r.manifest            = manifest.text;
    r.correlator_type     = {{"A_d-A_d", 'd'}, {"n_d-n_d", 's'}};
    r.nrg_revision        = "donrolih/nrgljubljana@c928cc55";
    r.z                   = {0.5, 1.0};
    r.complete_z_average  = true;
    for (double z : r.z) {
      twist_result_t t;
      t.z     = z;
      t.param = "[param]\nz=" + std::to_string(z) + "\n";
      t.spectra["A_d-A_d"] = nrg::spectrum_bins_t{{-0.1, 0.1}, {0.4, 0.6 * z}, {}, false};
      t.spectra["x-x"]     = nrg::spectrum_bins_t{{0.2}, {1.0}, {-0.5}, true};
      t.custom             = {{"T", "<n_d>"}, nda::array<double, 2>{{1.0, 0.9}, {0.5, 1.0 * z}}};
      t.td                 = {{"T", "E"}, nda::array<double, 2>{{1.0, -0.3}}};
      t.customfdm          = {{"n_d", z}};
      t.tdfdm              = {{"F", -z}};
      t.energy_flow        = {{{0.0, {{0, 1}}, 1}, {0.4, {{-1, 2}, {1, 2}}, 4}}, {}, {{0.1, {{0, 3}}, 3}}};
      r.twists.push_back(t);
    }
    r.expv  = {{"n_d", 0.75}};
    r.tdfdm = {{"F", -0.75}};
    auto mesh = triqs::mesh::refreq_log{1e-3, 1.0, 1.5};
    r.mesh    = mesh;
    g_w_t g{mesh, {{"imp", 1}}};
    for (auto w : mesh) g[0][w](0, 0) = 1.0 / (double(w) + 0.1i);
    r.processed["G"]     = g;
    r.processed["Sigma"] = g;
    return r;
  }

  void expect_same(solve_result_t const &a, solve_result_t const &b, bool with_bins) {
    EXPECT_EQ(a.request, b.request);
    EXPECT_EQ(a.solve_params.T, b.solve_params.T);
    EXPECT_EQ(a.solve_params.Nz, b.solve_params.Nz);
    EXPECT_EQ(a.solve_params.model_parameters, b.solve_params.model_parameters);
    EXPECT_EQ(a.nrg_params.bins, b.nrg_params.bins);
    EXPECT_EQ(a.model, b.model);
    EXPECT_EQ(a.symtype, b.symtype);
    EXPECT_EQ(a.manifest, b.manifest);
    EXPECT_EQ(a.correlator_type, b.correlator_type);
    EXPECT_EQ(a.nrg_revision, b.nrg_revision);
    EXPECT_EQ(a.mesh, b.mesh);
    EXPECT_EQ(a.z, b.z);
    EXPECT_EQ(a.complete_z_average, b.complete_z_average);
    ASSERT_EQ(a.twists.size(), b.twists.size());
    for (std::size_t k = 0; k < a.twists.size(); ++k) {
      auto const &s = a.twists[k], &t = b.twists[k];
      EXPECT_EQ(s.z, t.z);
      EXPECT_EQ(s.param, t.param);
      if (with_bins)
        EXPECT_EQ(s.spectra, t.spectra);
      else
        EXPECT_TRUE(t.spectra.empty());
      EXPECT_EQ(s.custom.columns, t.custom.columns);
      EXPECT_ARRAY_EQ(s.custom.values, t.custom.values);
      EXPECT_EQ(s.td.columns, t.td.columns);
      EXPECT_ARRAY_EQ(s.td.values, t.td.values);
      EXPECT_EQ(s.customfdm, t.customfdm);
      EXPECT_EQ(s.tdfdm, t.tdfdm);
      EXPECT_EQ(s.energy_flow, t.energy_flow);
    }
    EXPECT_EQ(a.expv, b.expv);
    EXPECT_EQ(a.tdfdm, b.tdfdm);
    ASSERT_EQ(a.processed.size(), b.processed.size());
    for (auto const &[name, g] : a.processed) EXPECT_BLOCK_GF_NEAR(g, b.processed.at(name), 0.0);
  }

} // namespace

TEST(SolveResult, H5RoundTrip) {
  auto const r = sample_result();
  {
    h5::file f("solve_result.h5", 'w');
    h5_write(f, "with", r);
    h5_write(f, "without", r, false);
  }
  h5::file f("solve_result.h5", 'r');
  expect_same(r, h5::read<solve_result_t>(f, "with"), true);
  expect_same(r, h5::read<solve_result_t>(f, "without"), false);
}

MAKE_MAIN
