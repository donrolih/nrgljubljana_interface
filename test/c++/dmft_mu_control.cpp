// Chemical potential (dmft/loop/mu_control) and hybridization_from_gamma (dmft/causality). No NRG run.
//
// The lattice is a tabulated smooth band, ρ ∝ (1 − ε²)² on [−1, 1], with Σ = −10⁻⁴ i: its spectral function is C¹ at
// the band edges, so the spectral weight on the mesh is 1 well within weight_tol. The filling n(μ) of the same
// lattice and occupancy, found by bisection, is the reference for the root.

#include <gtest/gtest.h>

#include <nrgljubljana_interface/dmft/causality.hpp>
#include <nrgljubljana_interface/dmft/lattice/tabulated_lattice.hpp>
#include <nrgljubljana_interface/dmft/loop/mu_control.hpp>
#include <nrgljubljana_interface/dmft/measurements.hpp>
#include <nrgljubljana_interface/nrg/integrate.hpp>

#include <cmath>
#include <complex>
#include <filesystem>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>

using namespace nrgljubljana_interface;
using namespace nrgljubljana_interface::dmft;
using triqs::mesh::refreq_log;
using dcomplex = std::complex<double>;

namespace {

  refreq_log test_mesh() { return {1e-4, 10.0, 1.05}; }

  dos_t smooth_band() {
    dos_t d;
    long const n = 201;
    for (long j = 0; j < n; ++j) {
      double const e = -1.0 + 2.0 * double(j) / double(n - 1);
      d.energy.push_back(e);
      d.density.push_back(std::pow(1 - e * e, 2));
    }
    double const norm = nrg::integrate(d.energy, d.density);
    for (auto &r : d.density) r /= norm;
    return d;
  }

  g_w_t sigma(refreq_log const &mesh) {
    g_w_t S{mesh, triqs::gfs::gf_struct_t{{"imp", 1}}};
    S[0].data() = dcomplex(0, -1e-4);
    return S;
  }

  double const T  = 0.01;
  double const mu0 = 0.1;
  std::map<std::string, long> const degeneracy{{"imp", 2}};

  struct setup_t {
    refreq_log mesh = test_mesh();
    tabulated_lattice_t lattice{smooth_band()};
    g_w_t Sigma = sigma(mesh);

    double filling(double mu) const { return occupancy(lattice.local_gf(Sigma, mu), T, degeneracy).total; }

    // The μ with n(μ) = goal, by bisection in [lo, hi]
    double bisection(double goal, double lo, double hi) const {
      for (int i = 0; i < 32; ++i) {
        double const mid = (lo + hi) / 2;
        (filling(mid) < goal ? lo : hi) = mid;
      }
      return (lo + hi) / 2;
    }

    mu_input_t input(double mu = mu0) const { return {Sigma, mu, T, degeneracy}; }
  };

  mu_control_params_t with_goal(double goal) {
    mu_control_params_t p;
    p.goal = goal;
    return p;
  }

  template <typename Exception, typename Function> void expect_error(Function f, std::string const &detail) {
    try {
      f();
    } catch (Exception const &e) {
      EXPECT_NE(std::string(e.what()).find(detail), std::string::npos) << "message: " << e.what();
      return;
    }
    ADD_FAILURE() << "no exception containing: " << detail;
  }

} // namespace

TEST(MuControl, Fixed) {
  setup_t const s;
  auto const adjusted = mu_control_t{}.adjust(s.input(), s.lattice);
  EXPECT_EQ(adjusted.mu_next, mu0);
  EXPECT_EQ(adjusted.report.status, "fixed");
  EXPECT_EQ(adjusted.filling, s.filling(mu0));
  EXPECT_EQ(adjusted.report.evaluations, 1);
  EXPECT_LT(adjusted.report.weight_error_used, 1e-4);
}

TEST(MuControl, WithinTolerance) {
  setup_t const s;
  auto const adjusted = mu_control_t{with_goal(s.filling(mu0))}.adjust(s.input(), s.lattice);
  EXPECT_EQ(adjusted.mu_next, mu0);
  EXPECT_EQ(adjusted.report.status, "within_tolerance");
  EXPECT_EQ(adjusted.report.evaluations, 1);
}

TEST(MuControl, Bracketed) {
  setup_t const s;
  for (double target : {mu0 + 0.02, mu0 - 0.03}) {
    double const goal   = s.filling(target);
    auto const adjusted = mu_control_t{with_goal(goal)}.adjust(s.input(), s.lattice);
    EXPECT_EQ(adjusted.report.status, "bracketed") << "target " << target;
    ASSERT_TRUE(adjusted.report.mu_root.has_value());
    double const reference = s.bisection(goal, mu0 - 0.1, mu0 + 0.1);
    EXPECT_NEAR(*adjusted.report.mu_root, reference, 1e-6) << "target " << target;
    EXPECT_NEAR(adjusted.mu_next, *adjusted.report.mu_root, 1e-15);
    EXPECT_LE(std::abs(adjusted.report.filling_next - goal), 1e-6);
    EXPECT_LE(adjusted.report.evaluations, 12);
  }
}

TEST(MuControl, Limited) {
  setup_t const s;
  for (double sign : {1.0, -1.0}) {
    auto const adjusted = mu_control_t{with_goal(s.filling(mu0 + sign * 0.3))}.adjust(s.input(), s.lattice);
    EXPECT_EQ(adjusted.report.status, "limited");
    EXPECT_NEAR(adjusted.mu_next, mu0 + sign * 0.05, 1e-15);
    EXPECT_LE(std::abs(adjusted.mu_next - mu0), 0.05);
    EXPECT_FALSE(adjusted.report.mu_root.has_value());
  }
}

// under = 0.5: the bracket reaches maxdx / under = 0.1, the step is half the way to the root
TEST(MuControl, Damped) {
  setup_t const s;
  auto p              = with_goal(s.filling(mu0 + 0.06));
  p.under             = 0.5;
  auto const adjusted = mu_control_t{p}.adjust(s.input(), s.lattice);
  EXPECT_EQ(adjusted.report.status, "bracketed");
  ASSERT_TRUE(adjusted.report.mu_root.has_value());
  EXPECT_NEAR(adjusted.mu_next, mu0 + 0.5 * (*adjusted.report.mu_root - mu0), 1e-15);
  EXPECT_NEAR(*adjusted.report.mu_root, mu0 + 0.06, 1e-5);
}

TEST(MuControl, Budget) {
  setup_t const s;
  auto p            = with_goal(s.filling(mu0 + 0.02));
  p.max_evaluations = 4;
  p.mu_tol          = 1e-15;
  expect_error<std::runtime_error>([&] { (void)mu_control_t{p}.adjust(s.input(), s.lattice); }, "max_evaluations (4)");
}

TEST(MuControl, WeightCheck) {
  setup_t const s;
  g_w_t narrow = sigma(refreq_log{1e-4, 0.6, 1.05}); // the mesh cuts the band
  mu_input_t in{narrow, mu0, T, degeneracy};
  expect_error<std::runtime_error>([&] { (void)mu_control_t{with_goal(0.8)}.adjust(in, s.lattice); }, "the lattice spectral weight at μ = 0.1 is off by");
  auto const fixed = mu_control_t{}.adjust(in, s.lattice); // without a goal the weight error is only reported
  EXPECT_GT(fixed.report.weight_error_used, 0.1);
}

TEST(MuControl, Parameters) {
  auto bad = [](auto change, std::string const &detail) {
    mu_control_params_t p;
    change(p);
    expect_error<std::invalid_argument>([&] { mu_control_t{p}; }, detail);
  };
  bad([](auto &p) { p.goal = -0.1; }, "the goal must be finite and >= 0");
  bad([](auto &p) { p.maxdx = 0; }, "maxdx must be positive");
  bad([](auto &p) { p.under = 1.5; }, "under must satisfy 0 < under ≤ 1");
  bad([](auto &p) { p.root_tol = 0; }, "root_tol must be positive");
  bad([](auto &p) { p.max_evaluations = 3; }, "max_evaluations must be at least 4");
  bad([](auto &p) { p.weight_tol = std::numeric_limits<double>::quiet_NaN(); }, "weight_tol must be positive");

  std::string const path = "dmft_mu_control.h5";
  auto p                 = with_goal(0.8);
  p.maxdx                = 0.02;
  p.max_evaluations      = 9;
  {
    h5::file f(path, 'w');
    h5_write(h5::group(f), "mu", p);
    h5_write(h5::group(f), "fixed", mu_control_params_t{});
  }
  mu_control_params_t r, q;
  {
    h5::file f(path, 'r');
    h5_read(h5::group(f), "mu", r);
    h5_read(h5::group(f), "fixed", q);
  }
  EXPECT_EQ(r.goal, p.goal);
  EXPECT_EQ(r.maxdx, p.maxdx);
  EXPECT_EQ(r.max_evaluations, p.max_evaluations);
  EXPECT_FALSE(q.goal.has_value());
  std::filesystem::remove(path);
}

// The step from a mixed Γ to Δ equals the projection of −iΓ followed by the rebuild
TEST(MuControl, HybridizationFromGamma) {
  refreq_log const mesh = test_mesh();
  g_w_t Gamma{mesh, triqs::gfs::gf_struct_t{{"imp", 1}, {"up", 2}}};
  auto const &pts = mesh.points();
  for (long k = 0; k < long(pts.size()); ++k) {
    double const w           = pts[k];
    Gamma[0].data()(k, 0, 0) = 0.3 / (1 + 4 * w * w) - 0.01; // negative at large |ω|
    Gamma[1].data()(k, 0, 0) = 0.5 / (1 + w * w);
    Gamma[1].data()(k, 1, 1) = 0.2 / (1 + (w - 0.5) * (w - 0.5));
    Gamma[1].data()(k, 0, 1) = dcomplex(0.25, 0.1) / (1 + w * w);
    Gamma[1].data()(k, 1, 0) = std::conj(Gamma[1].data()(k, 0, 1));
  }
  block_matrices_t h{nda::matrix<dcomplex>{{0.1}}, nda::matrix<dcomplex>{{-0.2, dcomplex(0.05, 0.01)}, {dcomplex(0.05, -0.01), 0.3}}};
  double const floor   = 1e-6;
  auto const [D, rep]  = hybridization_from_gamma(Gamma, h, floor);
  g_w_t minus_i_gamma  = Gamma;
  for (long b = 0; b < 2; ++b) minus_i_gamma[b].data() = dcomplex(0, -1) * Gamma[b].data();
  auto const [P, rep2] = causal_projection(minus_i_gamma, floor);
  auto const R         = rebuild_hermitian_part(P, h);
  EXPECT_EQ(rep, rep2);
  EXPECT_GT(rep.projected_points, 0);
  for (long b = 0; b < 2; ++b) EXPECT_LT(double(max_element(abs(D[b].data() - R[b].data()))), 1e-13) << "block " << b;

  Gamma[1].data()(5, 0, 1) += 1e-3;
  expect_error<std::invalid_argument>([&] { (void)hybridization_from_gamma(Gamma, h, floor); }, "Gamma is not Hermitian in block up");
}
