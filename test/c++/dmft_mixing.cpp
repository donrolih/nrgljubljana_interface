// Mixing (dmft/loop/mixing). No NRG run.

#include <gtest/gtest.h>

#include <nrgljubljana_interface/dmft/loop/mixing.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace nrgljubljana_interface;
using namespace nrgljubljana_interface::dmft;
using triqs::mesh::refreq_log;
using dcomplex = std::complex<double>;

namespace {

  refreq_log test_mesh() { return {1e-2, 5.0, 1.2}; }

  // Hermitian Γ in blocks imp (1x1) and up (2x2), varying with ω and with a parameter
  g_w_t make_Gamma(refreq_log const &mesh, double a) {
    g_w_t g{mesh, triqs::gfs::gf_struct_t{{"imp", 1}, {"up", 2}}};
    auto const &pts = mesh.points();
    for (long k = 0; k < long(pts.size()); ++k) {
      double const w       = pts[k];
      g[0].data()(k, 0, 0) = a / (1 + w * w);
      g[1].data()(k, 0, 0) = a * 0.5 / (1 + (w - 0.3) * (w - 0.3));
      g[1].data()(k, 1, 1) = 0.2 + a * 0.1 * std::exp(-w * w);
      g[1].data()(k, 0, 1) = dcomplex(0.05 * a, 0.02 * w);
      g[1].data()(k, 1, 0) = std::conj(g[1].data()(k, 0, 1));
    }
    return g;
  }

  nda::matrix<dcomplex> make_matrix_zero(long n) {
    nda::matrix<dcomplex> m(n, n);
    m = 0;
    return m;
  }

  mix_input_t input(refreq_log const &mesh) { return {3, make_Gamma(mesh, 1.0), make_Gamma(mesh, 2.5), 0.4, 0.8}; }

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

TEST(Mixing, LinearUpdate) {
  auto const mesh = test_mesh();
  auto const in   = input(mesh);
  double const a  = 0.3;
  linear_mixer mixer(a);
  auto const mixed = mixer.mix(in);
  EXPECT_FALSE(mixed.mu_next.has_value());
  EXPECT_EQ(mixed.report.method, "linear");
  for (long b = 0; b < 2; ++b) {
    auto const &out = mixed.Gamma_mix[b].data();
    auto const &raw = in.Gamma_out[b].data();
    auto const &old = in.Gamma[b].data();
    for (long k = 0; k < out.extent(0); ++k)
      for (long i = 0; i < out.extent(1); ++i)
        for (long j = 0; j < out.extent(2); ++j) {
          EXPECT_EQ(out(k, i, j), a * raw(k, i, j) + (1 - a) * old(k, i, j));
          EXPECT_EQ(out(k, i, j), std::conj(out(k, j, i)));
        }
  }
  // α = 1 gives Γ_out
  auto const full = linear_mixer(1.0).mix(in);
  for (long b = 0; b < 2; ++b) EXPECT_EQ(double(max_element(abs(full.Gamma_mix[b].data() - in.Gamma_out[b].data()))), 0.0);
}

TEST(Mixing, ResidualNorm) {
  auto const mesh = test_mesh();
  g_w_t a{mesh, triqs::gfs::gf_struct_t{{"imp", 1}}};
  g_w_t b = a;
  a[0].data()  = 0.7;
  b[0].data()  = 0.4;
  auto const &pts     = mesh.points();
  double const range  = pts.back() - pts.front();
  auto const mixed    = linear_mixer(0.5).mix({1, b, a, 0.0, 0.0});
  EXPECT_NEAR(mixed.report.residual_norm, 0.3 * std::sqrt(range), 1e-12);
}

// mix_hybridization is mix on the spectral parts followed by hybridization_from_Gamma
TEST(Mixing, MixHybridization) {
  auto const mesh  = test_mesh();
  block_matrices_t const onsite{make_matrix_zero(1), make_matrix_zero(2)};
  auto const Delta     = hybridization_from_Gamma(make_Gamma(mesh, 1.0), onsite, 1e-6).first;
  auto const Delta_out = hybridization_from_Gamma(make_Gamma(mesh, 2.5), onsite, 1e-6).first;
  linear_mixer a(0.3), b(0.3);
  auto const mixed = a.mix_hybridization(4, Delta, Delta_out, onsite, 1e-6);
  auto const by_hand = b.mix({4, spectral_part(Delta), spectral_part(Delta_out), 0.0, 0.0});
  auto const [Delta_mix, projection] = hybridization_from_Gamma(by_hand.Gamma_mix, onsite, 1e-6);
  for (long bl = 0; bl < 2; ++bl) EXPECT_EQ(double(max_element(abs(mixed.Delta_mix[bl].data() - Delta_mix[bl].data()))), 0.0);
  EXPECT_EQ(mixed.projection, projection);
  EXPECT_EQ(mixed.report.method, "linear");
  EXPECT_FALSE(mixed.mu_next.has_value());
}

TEST(Mixing, CheckpointRoundTrip) {
  auto const path = std::string("dmft_mixing.h5");
  {
    std::unique_ptr<mixer_t> m = std::make_unique<linear_mixer>(0.4);
    h5::file f(path, 'w');
    m->h5_write(h5::group(f), "mixer");
  }
  {
    std::unique_ptr<mixer_t> m = std::make_unique<linear_mixer>(0.4);
    h5::file f(path, 'r');
    m->h5_read(h5::group(f), "mixer");
    double alpha = 0;
    h5::read(h5::group(f).open_group("mixer"), "alpha", alpha);
    EXPECT_EQ(alpha, 0.4);
  }
  {
    h5::file f(path, 'w');
    auto g = h5::group(f).create_group("mixer");
    h5::write(g, "type", std::string("broyden"));
  }
  {
    linear_mixer m(0.4);
    h5::file f(path, 'r');
    expect_error<std::runtime_error>([&] { m.h5_read(h5::group(f), "mixer"); }, "linear mixer: the group mixer holds a broyden mixer");
  }
  std::filesystem::remove(path);
}

TEST(Mixing, Errors) {
  auto const mesh = test_mesh();
  expect_error<std::invalid_argument>([] { linear_mixer(0.0); }, "α must satisfy 0 < α ≤ 1, got 0");
  expect_error<std::invalid_argument>([] { linear_mixer(1.5); }, "got 1.5");
  expect_error<std::invalid_argument>([] { linear_mixer(std::numeric_limits<double>::quiet_NaN()); }, "α must satisfy");

  linear_mixer mixer(0.5);
  auto in = input(mesh);
  g_w_t other{mesh, triqs::gfs::gf_struct_t{{"imp", 1}}};
  expect_error<std::invalid_argument>([&] { mixer.mix({1, in.Gamma, other, 0.0, 0.0}); }, "Gamma_out does not have the blocks of Gamma");
  g_w_t other_mesh = make_Gamma(refreq_log{1e-2, 5.0, 1.3}, 1.0);
  expect_error<std::invalid_argument>([&] { mixer.mix({1, in.Gamma, other_mesh, 0.0, 0.0}); }, "is not on the mesh of Gamma");
  in.Gamma_out[1].data()(2, 1, 1) = std::numeric_limits<double>::infinity();
  expect_error<std::invalid_argument>([&] { mixer.mix(in); }, "Gamma_out is not finite in block up, element (1, 1)");
}

// ---------------------------------------------------------------- Broyden

namespace {

  // A positive definite Γ in blocks imp (1x1) and up (2x2), varying with ω and with a parameter
  g_w_t positive_Gamma(refreq_log const &mesh, double a) {
    g_w_t g{mesh, triqs::gfs::gf_struct_t{{"imp", 1}, {"up", 2}}};
    auto const &pts = mesh.points();
    for (long k = 0; k < long(pts.size()); ++k) {
      double const w       = pts[k];
      g[0].data()(k, 0, 0) = 0.1 + a / (1 + w * w);
      g[1].data()(k, 0, 0) = 0.3 + a * 0.5 / (1 + (w - 0.3) * (w - 0.3));
      g[1].data()(k, 1, 1) = 0.2 + a * 0.1 * std::exp(-w * w);
      g[1].data()(k, 0, 1) = dcomplex(0.02 * a, 0.01 * w / (1 + w * w));
      g[1].data()(k, 1, 0) = std::conj(g[1].data()(k, 0, 1));
    }
    return g;
  }

  // A linear fixed-point problem: Γ_out = Γ* + m (Γ − Γ*), with a factor m per block and half of the mesh (four
  // distinct factors), so that a multisecant method solves it in a few steps
  struct linear_problem {
    refreq_log mesh = test_mesh();
    g_w_t star      = positive_Gamma(mesh, 2.0);

    double factor(long b, double w) const { return (b == 0 ? 0.8 : 0.3) * (w < 0 ? 1.0 : -0.5); }

    g_w_t out(g_w_t const &Gamma) const {
      g_w_t r         = Gamma;
      auto const &pts = mesh.points();
      for (long b = 0; b < 2; ++b)
        for (long k = 0; k < long(pts.size()); ++k) {
          double const m = factor(b, pts[k]);
          for (long i = 0; i < r[b].data().extent(1); ++i)
            for (long j = 0; j < r[b].data().extent(2); ++j)
              r[b].data()(k, i, j) = star[b].data()(k, i, j) + m * (Gamma[b].data()(k, i, j) - star[b].data()(k, i, j));
        }
      return r;
    }

    // The residual norms of the iterations, until the norm is below reduction times the first one (at most
    // max_iterations)
    std::vector<double> run(mixer_t &mixer, long max_iterations, double reduction) const {
      g_w_t Gamma = positive_Gamma(mesh, 2.3);
      std::vector<double> norms;
      for (long it = 1; it <= max_iterations; ++it) {
        auto const mixed = mixer.mix({it, Gamma, out(Gamma), 0.0, 0.0});
        norms.push_back(mixed.report.residual_norm);
        if (mixed.report.residual_norm < reduction * norms.front()) break;
        Gamma = mixed.Gamma_mix;
      }
      return norms;
    }
  };

  // A constant scalar Γ in one block
  g_w_t constant_Gamma(double value) {
    g_w_t g{test_mesh(), triqs::gfs::gf_struct_t{{"imp", 1}}};
    g[0].data() = value;
    return g;
  }

} // namespace

TEST(Broyden, LinearFixedPoint) {
  linear_problem const problem;
  broyden_mixer broyden;
  auto const norms = problem.run(broyden, 15, 1e-6);
  EXPECT_LT(norms.back(), 1e-6 * norms.front());
  EXPECT_EQ(broyden.history_size(), std::min<long>(10, long(norms.size())));

  // linear mixing reduces the residual by at most |1 − α(1 − 0.8)| = 0.9 per step
  linear_mixer linear(0.5);
  auto const linear_norms = problem.run(linear, long(norms.size()), 0);
  EXPECT_GT(linear_norms.back(), 1e-3 * linear_norms.front());

  // unweighted, as dmft-jamova
  broyden_params_t p;
  p.weighted = false;
  broyden_mixer unweighted(p);
  auto const u = problem.run(unweighted, 15, 1e-6);
  EXPECT_LT(u.back(), 1e-6 * u.front());
}

TEST(Broyden, Mu) {
  // Γ at its fixed point; n(μ) = 0.8 + 0.7 (μ − 0.53): the secant finds the root after the first step
  auto const Gamma = constant_Gamma(0.4);
  auto filling     = [](double mu) { return 0.8 + 0.7 * (mu - 0.53); };
  broyden_params_t p;
  p.mu     = true;
  p.warmup = 2;
  broyden_mixer mixer(p);
  EXPECT_FALSE(mixer.moves_mu(2));
  EXPECT_TRUE(mixer.moves_mu(3));
  double mu = 0.5;
  for (long it = 1; it <= 2; ++it) {
    auto const mixed = mixer.mix({it, Gamma, Gamma, mu, filling(mu), 0.8});
    EXPECT_FALSE(mixed.mu_next.has_value());
    mu += 0.01; // as a μ control would move it
  }
  for (long it = 3; it <= 5; ++it) {
    auto const mixed = mixer.mix({it, Gamma, Gamma, mu, filling(mu), 0.8});
    ASSERT_TRUE(mixed.mu_next.has_value());
    mu = *mixed.mu_next;
  }
  EXPECT_NEAR(mu, 0.53, 1e-6);
}

TEST(Broyden, MuStepClamped) {
  auto const Gamma = constant_Gamma(0.4);
  broyden_params_t p;
  p.mu          = true;
  p.mu_alpha    = 1.0;
  p.max_mu_step = 0.05;
  broyden_mixer mixer(p);
  auto const mixed = mixer.mix({1, Gamma, Gamma, 0.0, 0.8 + 0.7 * (0.0 - 0.53), 0.8}); // the linear step is 0.371
  ASSERT_TRUE(mixed.mu_next.has_value());
  EXPECT_DOUBLE_EQ(*mixed.mu_next, 0.05);
}

TEST(Broyden, Fallbacks) {
  // a repeated residual gives no secant pair: the linear step, not counted as a fallback
  {
    broyden_mixer mixer;
    auto const G = constant_Gamma(1.0), G_out = constant_Gamma(1.2);
    (void)mixer.mix({1, G, G_out, 0, 0});
    auto const mixed = mixer.mix({2, G, G_out, 0, 0});
    EXPECT_FALSE(mixed.report.fallback);
    EXPECT_NEAR(double(mixed.Gamma_mix[0].data()(0, 0, 0).real()), 1.1, 1e-14);
  }
  // an ill-conditioned system
  {
    broyden_params_t p;
    p.max_condition = 1.0;
    broyden_mixer mixer(p);
    linear_problem const problem;
    g_w_t G = positive_Gamma(problem.mesh, 2.3);
    mix_result_t mixed;
    for (long it = 1; it <= 3; ++it) {
      mixed = mixer.mix({it, G, problem.out(G), 0, 0});
      G     = mixed.Gamma_mix;
    }
    EXPECT_TRUE(mixed.report.fallback);
    EXPECT_EQ(mixed.report.method, "linear (fallback)");
  }
  // a step to a negative Γ: secant root at −0.25, linear step 0.67
  {
    broyden_mixer mixer;
    (void)mixer.mix({1, constant_Gamma(1.0), constant_Gamma(0.5), 0, 0});
    auto const mixed = mixer.mix({2, constant_Gamma(0.9), constant_Gamma(0.44), 0, 0});
    EXPECT_TRUE(mixed.report.fallback);
    EXPECT_NEAR(double(mixed.Gamma_mix[0].data()(0, 0, 0).real()), 0.67, 1e-12);
  }
}

// Interrupted after three iterations, written and read into a new mixer: the same steps as without
TEST(Broyden, Checkpoint) {
  linear_problem const problem;
  std::string const path = "dmft_mixing_broyden.h5";
  broyden_mixer full;
  g_w_t G = positive_Gamma(problem.mesh, 2.3);
  std::vector<g_w_t> steps;
  for (long it = 1; it <= 6; ++it) {
    G = full.mix({it, G, problem.out(G), 0, 0}).Gamma_mix;
    steps.push_back(G);
  }

  broyden_mixer first;
  G = positive_Gamma(problem.mesh, 2.3);
  for (long it = 1; it <= 3; ++it) G = first.mix({it, G, problem.out(G), 0, 0}).Gamma_mix;
  {
    h5::file f(path, 'w');
    first.h5_write(h5::group(f), "mixer");
  }
  broyden_mixer second;
  {
    h5::file f(path, 'r');
    second.h5_read(h5::group(f), "mixer");
  }
  EXPECT_EQ(second.history_size(), 3);
  for (long it = 4; it <= 6; ++it) {
    G = second.mix({it, G, problem.out(G), 0, 0}).Gamma_mix;
    for (long b = 0; b < 2; ++b) EXPECT_EQ(double(max_element(abs(G[b].data() - steps[it - 1][b].data()))), 0.0);
  }
  {
    h5::file f(path, 'r');
    linear_mixer linear(0.5);
    expect_error<std::runtime_error>([&] { linear.h5_read(h5::group(f), "mixer"); }, "holds a broyden mixer");
  }
  std::filesystem::remove(path);
}

TEST(Broyden, Errors) {
  auto bad = [](auto change, std::string const &detail) {
    broyden_params_t p;
    change(p);
    expect_error<std::invalid_argument>([&] { broyden_mixer{p}; }, detail);
  };
  bad([](auto &p) { p.alpha = 0; }, "alpha must satisfy 0 < alpha ≤ 1");
  bad([](auto &p) { p.history = 0; }, "history must be at least 1");
  bad([](auto &p) { p.mu_alpha = -1; }, "mu_alpha must be positive");
  bad([](auto &p) { p.warmup = -1; }, "warmup must be >= 0");
  bad([](auto &p) { p.max_mu_step = 0; }, "max_mu_step must be positive");
  bad([](auto &p) { p.max_condition = 0.5; }, "max_condition must be at least 1");

  auto const G = constant_Gamma(1.0);
  broyden_mixer mixer;
  (void)mixer.mix({1, G, G, 0, 0});
  expect_error<std::invalid_argument>([&] { (void)mixer.mix({3, G, G, 0, 0}); }, "iteration 3 does not follow iteration 1");
  broyden_params_t p;
  p.mu = true;
  broyden_mixer with_mu(p);
  expect_error<std::invalid_argument>([&] { (void)with_mu.mix({1, G, G, 0, 0}); }, "mixing μ needs the goal of the filling");
  EXPECT_NE(broyden_mixer{}.description().find("broyden alpha=0.5 history=10"), std::string::npos);
}
