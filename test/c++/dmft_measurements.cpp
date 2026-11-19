// Measurements (dmft/measurements). No NRG run.
//
// The spectra are compact bands, A(ω) = S_B(ω)/π with S_B(ω) = G0 (1 − ((ω − c)/D)²)² for |ω − c| < D, whose
// integrals have closed forms. On the mesh used here the Steffen interpolant integrates the band weight to about 1e-6
// relative (estimated with the same interpolant in numpy); the tolerances against closed forms are 1e-5. Matrix cases
// are compared with the scalar result times the matrix: the interpolant is linear in the data.

#include <gtest/gtest.h>

#include <nrgljubljana_interface/dmft/measurements.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <map>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

using namespace nrgljubljana_interface;
using namespace nrgljubljana_interface::dmft;
using triqs::mesh::refreq_log;
using dcomplex = std::complex<double>;
using cmatrix  = nda::matrix<dcomplex>;

namespace {

  refreq_log test_mesh() { return {1e-4, 10.0, 1.03}; }

  struct band_t {
    double c, D, G0;
  };

  // H_B − i S_B, with the closed form of H_B = KK(−S_B)
  dcomplex band(double w, band_t const &b) {
    double const u   = w - b.c;
    double const P   = b.G0 * std::pow(1 - (u / b.D) * (u / b.D), 2);
    double const lg  = std::abs(std::abs(u) - b.D) == 0 ? 0.0 : std::log(std::abs((u + b.D) / (u - b.D)));
    double const H   = (P * lg + b.G0 * (10 * u / (3 * b.D) - 2 * u * u * u / (b.D * b.D * b.D))) / std::numbers::pi;
    double const S_b = std::abs(u) < b.D ? P : 0.0;
    return {H, -S_b};
  }

  // ∫ A dω = ∫ S_B dω / π
  double band_weight(band_t const &b) { return 16 * b.D * b.G0 / (15 * std::numbers::pi); }

  // ∫_{−∞}^{0} A dω
  double band_filling_t0(band_t const &b) {
    auto F            = [&](double t) { return b.G0 * (t - 2 * t * t * t / (3 * b.D * b.D) + std::pow(t, 5) / (5 * std::pow(b.D, 4))); };
    double const lo   = -b.D, hi = std::min(-b.c, b.D);
    return hi > lo ? (F(hi) - F(lo)) / std::numbers::pi : 0.0;
  }

  cmatrix make_matrix(std::initializer_list<std::initializer_list<dcomplex>> rows) {
    long const n = long(rows.size());
    cmatrix m(n, n);
    long i = 0;
    for (auto const &r : rows) {
      long j = 0;
      for (auto const &v : r) m(i, j++) = v;
      ++i;
    }
    return m;
  }

  cmatrix M_up() { return make_matrix({{1.0, {0.2, 0.3}}, {{0.2, -0.3}, 0.5}}); }

  // Blocks imp = B1 (1x1) and up = M_up B2 (2x2)
  g_w_t make_gf(refreq_log const &mesh, band_t const &b1, band_t const &b2) {
    g_w_t g{mesh, triqs::gfs::gf_struct_t{{"imp", 1}, {"up", 2}}};
    auto const &pts = mesh.points();
    cmatrix const M = M_up();
    for (long k = 0; k < long(pts.size()); ++k) {
      g[0].data()(k, 0, 0) = band(pts[k], b1);
      dcomplex const v     = band(pts[k], b2);
      for (long i = 0; i < 2; ++i)
        for (long j = 0; j < 2; ++j) g[1].data()(k, i, j) = M(i, j) * v;
    }
    return g;
  }

  // One 1x1 block holding a band
  g_w_t scalar_gf(refreq_log const &mesh, band_t const &b) {
    g_w_t g{mesh, triqs::gfs::gf_struct_t{{"imp", 1}}};
    auto const &pts = mesh.points();
    for (long k = 0; k < long(pts.size()); ++k) g[0].data()(k, 0, 0) = band(pts[k], b);
    return g;
  }

  band_t const b1{0.3, 1.0, 0.7}, b2{0.1, 0.8, 0.5};

  template <typename Function> void expect_invalid_argument(Function f, std::string const &detail) {
    try {
      f();
    } catch (std::invalid_argument const &e) {
      EXPECT_NE(std::string(e.what()).find(detail), std::string::npos) << "message: " << e.what();
      return;
    }
    ADD_FAILURE() << "no std::invalid_argument containing: " << detail;
  }

} // namespace

TEST(Measurements, SpectralWeight) {
  auto const mesh = test_mesh();
  auto const w    = spectral_weight(make_gf(mesh, b1, b2));
  ASSERT_EQ(w.size(), 2u);
  EXPECT_NEAR(w[0](0, 0).real(), band_weight(b1), 1e-5 * band_weight(b1));
  EXPECT_EQ(w[0](0, 0).imag(), 0.0);
  double const w2 = spectral_weight(scalar_gf(mesh, b2))[0](0, 0).real();
  EXPECT_NEAR(w2, band_weight(b2), 1e-5 * band_weight(b2));
  for (long i = 0; i < 2; ++i)
    for (long j = 0; j < 2; ++j) EXPECT_LT(std::abs(w[1](i, j) - M_up()(i, j) * w2), 1e-12) << i << j;
}

TEST(Measurements, OccupancyAtZeroTemperature) {
  auto const mesh = test_mesh();
  auto const n    = occupancy(scalar_gf(mesh, b1), 0.0, {{"imp", 1}});
  EXPECT_NEAR(n.matrix[0](0, 0).real(), band_filling_t0(b1), 1e-5 * band_weight(b1));
  EXPECT_EQ(n.total, n.matrix[0](0, 0).real());
}

// A band symmetric about ω = 0 is half filled at any T, since f(ω) + f(−ω) = 1, up to the quadrature tolerance
TEST(Measurements, OccupancySymmetricBand) {
  auto const mesh      = test_mesh();
  band_t const sym{0.0, 1.0, 0.7};
  auto const G         = scalar_gf(mesh, sym);
  double const weight  = spectral_weight(G)[0](0, 0).real();
  for (double T : {1e-3, 0.05, 0.5}) EXPECT_NEAR(occupancy(G, T, {{"imp", 1}}).matrix[0](0, 0).real(), weight / 2, 1e-8) << "T = " << T;
}

// The matrix block is M times the scalar result, up to the absolute quadrature tolerance of integrate_fermi
TEST(Measurements, OccupancyMatrixAndTotal) {
  auto const mesh = test_mesh();
  double const T  = 0.02;
  auto const n    = occupancy(make_gf(mesh, b1, b2), T, {{"imp", 2}, {"up", 1}});
  double const n1 = occupancy(scalar_gf(mesh, b1), T, {{"imp", 1}}).matrix[0](0, 0).real();
  double const n2 = occupancy(scalar_gf(mesh, b2), T, {{"imp", 1}}).matrix[0](0, 0).real();
  EXPECT_NEAR(n.matrix[0](0, 0).real(), n1, 1e-9);
  for (long i = 0; i < 2; ++i)
    for (long j = 0; j < 2; ++j) EXPECT_LT(std::abs(n.matrix[1](i, j) - M_up()(i, j) * n2), 1e-9) << i << j;
  EXPECT_EQ(n.matrix[1](0, 1), std::conj(n.matrix[1](1, 0)));
  EXPECT_NEAR(n.total, 2 * n1 + (M_up()(0, 0).real() + M_up()(1, 1).real()) * n2, 1e-8);
}

// Bands differing only in amplitude: ∫ |ΔA| = |ΔG0| ∫ S_B / (π G0); off-diagonal elements add the real and imaginary parts
TEST(Measurements, SpectralDistance) {
  auto const mesh = test_mesh();
  band_t const a1{0.3, 1.0, 0.7}, a2{0.3, 1.0, 0.5}, c1{0.1, 0.8, 0.5}, c2{0.1, 0.8, 0.45};
  auto const d    = spectral_distance(make_gf(mesh, a1, c1), make_gf(mesh, a2, c2));
  double const da = band_weight({0.3, 1.0, 0.2});
  double const dc = band_weight({0.1, 0.8, 0.05});
  ASSERT_EQ(d.elements.size(), 2u);
  EXPECT_NEAR(d.elements[0](0, 0), da, 1e-5 * da);
  double largest = d.elements[0](0, 0);
  for (long i = 0; i < 2; ++i)
    for (long j = 0; j < 2; ++j) {
      double const expected = (std::abs(M_up()(i, j).real()) + std::abs(M_up()(i, j).imag())) * dc;
      EXPECT_NEAR(d.elements[1](i, j), expected, 1e-5 * dc) << i << j;
      largest = std::max(largest, d.elements[1](i, j));
    }
  EXPECT_EQ(d.max, largest);
  EXPECT_EQ(spectral_distance(make_gf(mesh, a1, c1), make_gf(mesh, a1, c1)).max, 0.0);
}

TEST(Measurements, Errors) {
  auto const mesh = test_mesh();
  auto const G    = make_gf(mesh, b1, b2);
  std::map<std::string, long> const deg{{"imp", 2}, {"up", 1}};

  expect_invalid_argument([&] { occupancy(G, -1e-3, deg); }, "temperature must be finite and >= 0");
  expect_invalid_argument([&] { occupancy(G, std::numeric_limits<double>::quiet_NaN(), deg); }, "temperature must be finite and >= 0");
  expect_invalid_argument([&] { occupancy(G, 0.01, {{"imp", 2}}); }, "no degeneracy for block up");
  expect_invalid_argument([&] { occupancy(G, 0.01, {{"imp", 2}, {"up", 1}, {"dn", 1}}); }, "degeneracy for block dn, not in G");
  expect_invalid_argument([&] { occupancy(G, 0.01, {{"imp", 0}, {"up", 1}}); }, "the degeneracy of block imp is 0");

  auto bad                  = G;
  bad[1].data()(5, 1, 1) = {0.0, std::numeric_limits<double>::quiet_NaN()};
  expect_invalid_argument([&] { occupancy(bad, 0.01, deg); }, "occupancy: G is not finite in block up, element (1, 1)");
  expect_invalid_argument([&] { spectral_weight(bad); }, "spectral_weight: G is not finite in block up");
  expect_invalid_argument([&] { spectral_distance(G, bad); }, "spectral_distance: b is not finite in block up");
  expect_invalid_argument([&] { spectral_distance(G, scalar_gf(mesh, b1)); }, "spectral_distance: b does not have the blocks of a");
}
