// Mixing (dmft/loop/mixing). No NRG run.

#include <gtest/gtest.h>

#include <nrgljubljana_interface/dmft/loop/mixing.hpp>

#include <cmath>
#include <complex>
#include <filesystem>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

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
