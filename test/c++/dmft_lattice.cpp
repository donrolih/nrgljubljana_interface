// Lattices (dmft/lattice). No NRG run.
//
// The Bethe closed forms are compared with direct integration of ρ(ε) εⁿ (z − ε)⁻¹, with ε = D cos θ, which makes the
// integrand smooth and periodic, so the trapezoid rule converges exponentially (Im z ≥ 0.05 here). The tabulated
// semicircle is sampled uniformly in θ, which resolves the band edges, as dmft-jamova's DOS.dat; with NRG Ljubljana's
// hilb tool (analytic, Steffen) that table agrees with the closed form to about 1e-9 away from the edges, so the
// tolerance against bethe_lattice_t is 1e-7.

#include <gtest/gtest.h>

#include <nrgljubljana_interface/dmft/lattice/bethe_lattice.hpp>
#include <nrgljubljana_interface/dmft/lattice/tabulated_lattice.hpp>
#include <nrgljubljana_interface/nrg/integrate.hpp>

#include <nda/linalg.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <limits>
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

  refreq_log test_mesh() { return {1e-2, 4.0, 1.15}; }

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

  cmatrix at(nda::array<dcomplex, 3> const &d, long k) {
    cmatrix m(d.extent(1), d.extent(2));
    for (long i = 0; i < d.extent(1); ++i)
      for (long j = 0; j < d.extent(2); ++j) m(i, j) = d(k, i, j);
    return m;
  }

  cmatrix id(long n) {
    cmatrix m(n, n);
    m = 0;
    for (long i = 0; i < n; ++i) m(i, i) = 1;
    return m;
  }

  // Causal self-energies: imp 1x1, up 2x2 with off-diagonal elements (Hermitian A, positive M), Im parts ≤ −0.05
  cmatrix sigma_imp(double w) { return make_matrix({{0.1 + 0.4 / (w - 0.5 + dcomplex(0, 0.6)) - dcomplex(0, 0.05)}}); }
  cmatrix sigma_up(double w) {
    cmatrix const A = make_matrix({{0.2, {0.05, 0.1}}, {{0.05, -0.1}, -0.1}});
    cmatrix const M = make_matrix({{0.3, {0.1, -0.05}}, {{0.1, 0.05}, 0.2}});
    return A + M * (1.0 / (w + 0.3 + dcomplex(0, 0.4))) - dcomplex(0, 0.05) * id(2);
  }
  cmatrix sigma_up_diagonal(double w) { return make_matrix({{sigma_imp(w)(0, 0), 0.0}, {0.0, 0.3 - dcomplex(0, 0.1) + 0.2 / (w + dcomplex(0, 0.5))}}); }

  template <typename F> g_w_t make_sigma(refreq_log const &mesh, F up) {
    g_w_t S{mesh, triqs::gfs::gf_struct_t{{"imp", 1}, {"up", 2}}};
    auto const &pts = mesh.points();
    for (long k = 0; k < long(pts.size()); ++k) {
      S[0].data()(k, 0, 0) = sigma_imp(pts[k])(0, 0);
      cmatrix const u      = up(pts[k]);
      for (long i = 0; i < 2; ++i)
        for (long j = 0; j < 2; ++j) S[1].data()(k, i, j) = u(i, j);
    }
    return S;
  }

  // G_n = ∫ ρ(ε) εⁿ (z − ε)⁻¹ dε for the semicircle of half-bandwidth D, z = (ω + μ) − Σ
  cmatrix bethe_quadrature(cmatrix const &z, double D, int moment) {
    long const n = z.extent(0), N = 2000;
    cmatrix sum(n, n);
    sum = 0;
    for (long j = 0; j < N; ++j) {
      double const theta = 2 * std::numbers::pi * double(j) / double(N);
      double const e     = D * std::cos(theta);
      double const w     = std::pow(std::sin(theta), 2) * (moment == 0 ? 1.0 : e);
      sum += w * nda::linalg::inv(cmatrix(z - e * id(n)));
    }
    return sum * (2.0 / std::numbers::pi) * (std::numbers::pi / double(N)); // ∫_0^π = half of the periodic sum
  }

  // The semicircle sampled uniformly in θ (ε = −D cos θ) and normalized, optionally shifted
  dos_t semicircle(double D, long points = 2601, double shift = 0) {
    dos_t d;
    for (long j = 0; j < points; ++j) {
      double const e = -D * std::cos(std::numbers::pi * double(j) / double(points - 1));
      d.energy.push_back(e + shift);
      d.density.push_back(j == 0 or j == points - 1 ? 0.0 : 2 / (std::numbers::pi * D * D) * std::sqrt(std::max(0.0, D * D - e * e)));
    }
    double const norm = nrg::integrate(d.energy, d.density);
    for (auto &r : d.density) r /= norm;
    return d;
  }

  double max_abs_diff(g_w_t const &a, g_w_t const &b) {
    double m = 0;
    for (long bl = 0; bl < long(a.size()); ++bl) m = std::max(m, double(max_element(abs(a[bl].data() - b[bl].data()))));
    return m;
  }

  template <typename Function> void expect_invalid_argument(Function f, std::string const &detail) {
    try {
      f();
    } catch (std::invalid_argument const &e) {
      EXPECT_NE(std::string(e.what()).find(detail), std::string::npos) << "message: " << e.what();
      return;
    }
    ADD_FAILURE() << "no std::invalid_argument containing: " << detail;
  }

  double const D  = 1.3;
  double const mu = 0.2;

} // namespace

TEST(Lattice, BetheAgainstQuadrature) {
  auto const mesh  = test_mesh();
  auto const Sigma = make_sigma(mesh, sigma_up);
  auto const [G0, G1] = bethe_lattice_t(D).moments(Sigma, mu);
  auto const &pts     = mesh.points();
  for (long b = 0; b < 2; ++b)
    for (long k = 0; k < long(pts.size()); ++k) {
      long const n   = b == 0 ? 1 : 2;
      cmatrix const z = (pts[k] + mu) * id(n) - at(Sigma[b].data(), k);
      EXPECT_LT(double(max_element(abs(at(G0[b].data(), k) - bethe_quadrature(z, D, 0)))), 1e-10) << "block " << b << ", k = " << k;
      EXPECT_LT(double(max_element(abs(at(G1[b].data(), k) - bethe_quadrature(z, D, 1)))), 1e-10) << "block " << b << ", k = " << k;
    }
}

TEST(Lattice, BetheSelfConsistency) {
  auto const mesh       = test_mesh();
  auto const Sigma      = make_sigma(mesh, sigma_up);
  bethe_lattice_t const lattice(D);
  auto const [G_loc, Delta] = lattice.self_consistency(Sigma, mu);
  EXPECT_EQ(max_abs_diff(G_loc, lattice.local_gf(Sigma, mu)), 0.0);
  EXPECT_LT(max_abs_diff(Delta, g_w_t(D * D / 4 * G_loc)), 1e-12);
  // G1 = z G0 − 1
  auto const [G0, G1] = lattice.moments(Sigma, mu);
  auto const &pts     = mesh.points();
  for (long k = 0; k < long(pts.size()); ++k) {
    cmatrix const z = (pts[k] + mu) * id(2) - at(Sigma[1].data(), k);
    EXPECT_LT(double(max_element(abs(at(G1[1].data(), k) - (z * at(G0[1].data(), k) - id(2))))), 1e-12);
  }
  auto const h = lattice.onsite_energy({{"imp", 1}, {"up", 2}});
  ASSERT_EQ(h.size(), 2u);
  EXPECT_EQ(double(max_element(abs(h[1]))), 0.0);
}

TEST(Lattice, TabulatedSemicircle) {
  auto const mesh = test_mesh();
  tabulated_lattice_t const dos(semicircle(D));
  bethe_lattice_t const bethe(D);
  for (auto const &Sigma : {make_sigma(mesh, sigma_up), make_sigma(mesh, sigma_up_diagonal)}) {
    auto const [G0, G1] = dos.moments(Sigma, mu);
    auto const [B0, B1] = bethe.moments(Sigma, mu);
    EXPECT_LT(max_abs_diff(G0, B0), 1e-7);
    EXPECT_LT(max_abs_diff(G1, B1), 1e-7);
    EXPECT_EQ(max_abs_diff(dos.local_gf(Sigma, mu), G0), 0.0);
  }
  EXPECT_LT(std::abs(dos.onsite_energy({{"imp", 1}})[0](0, 0)), 1e-12);
}

TEST(Lattice, OnsiteEnergyOfAShiftedDos) {
  tabulated_lattice_t const dos(semicircle(D, 2601, 0.3));
  auto const h = dos.onsite_energy({{"imp", 1}, {"up", 2}});
  EXPECT_NEAR(h[0](0, 0).real(), 0.3, 1e-10);
  EXPECT_NEAR(h[1](1, 1).real(), 0.3, 1e-10);
  EXPECT_EQ(h[1](0, 1), dcomplex(0));
}

TEST(Lattice, DosPerElement) {
  auto const mesh = test_mesh();
  double const D2 = 0.6;
  tabulated_lattice_t const dos(std::map<std::string, std::vector<dos_t>>{{"imp", {semicircle(D)}}, {"up", {semicircle(D), semicircle(D2, 2601, 0.1)}}});
  auto const Sigma    = make_sigma(mesh, sigma_up_diagonal);
  auto const [G0, G1] = dos.moments(Sigma, mu);
  // Orbital 1 of up: the semicircle of D2 centred at 0.1, i.e. Bethe with μ − 0.1
  auto const [B0, B1] = bethe_lattice_t(D2).moments(Sigma, mu - 0.1);
  auto const &pts     = mesh.points();
  for (long k = 0; k < long(pts.size()); ++k) EXPECT_LT(std::abs(G0[1].data()(k, 1, 1) - B0[1].data()(k, 1, 1)), 1e-7) << "k = " << k;
  auto const [C0, C1] = bethe_lattice_t(D).moments(Sigma, mu);
  for (long k = 0; k < long(pts.size()); ++k) EXPECT_LT(std::abs(G0[1].data()(k, 0, 0) - C0[1].data()(k, 0, 0)), 1e-7) << "k = " << k;
  auto const h = dos.onsite_energy({{"imp", 1}, {"up", 2}});
  EXPECT_NEAR(h[1](1, 1).real(), 0.1, 1e-10);
  EXPECT_NEAR(h[1](0, 0).real(), 0.0, 1e-12);

  expect_invalid_argument([&] { dos.moments(make_sigma(mesh, sigma_up), mu); }, "off-diagonal elements in block up");
  expect_invalid_argument([&] { dos.onsite_energy({{"imp", 1}, {"up", 3}}); }, "block up has 3 orbitals, the DOS map 2");
  expect_invalid_argument([&] { dos.onsite_energy({{"imp", 1}}); }, "a DOS for block up, not in the structure");
  expect_invalid_argument([&] { dos.onsite_energy({{"imp", 1}, {"up", 2}, {"dn", 2}}); }, "no DOS for block dn");
}

TEST(Lattice, ReadDos) {
  auto const d    = semicircle(D, 101);
  std::string const path = "dmft_lattice_dos.dat";
  {
    std::ofstream out(path);
    out.precision(17);
    out << "# energy density\n\n";
    for (std::size_t k = 0; k < d.energy.size(); ++k) out << d.energy[k] << " " << d.density[k] << "\n";
  }
  auto const r = tabulated_lattice_t::read_dos(path);
  EXPECT_EQ(r.energy, d.energy);
  EXPECT_EQ(r.density, d.density);
  {
    std::ofstream out(path, std::ios::app);
    out << "1.5 0 7\n";
  }
  expect_invalid_argument([&] { tabulated_lattice_t::read_dos(path); }, "line 104: expected two numbers");
  std::filesystem::remove(path);
}

TEST(Lattice, Errors) {
  auto const mesh = test_mesh();
  expect_invalid_argument([] { bethe_lattice_t(0.0); }, "half-bandwidth must be positive and finite");
  expect_invalid_argument([] { bethe_lattice_t(1.0, 0.0); }, "the floor of Im Σ must be finite and at least");

  auto unnormalized = semicircle(D);
  for (auto &r : unnormalized.density) r *= 1.01;
  expect_invalid_argument([&] { tabulated_lattice_t{unnormalized}; }, "is not normalized: its integral is 1.01");
  auto negative          = semicircle(D);
  negative.density[1000] = -1e-3;
  expect_invalid_argument([&] { tabulated_lattice_t{negative}; }, "is negative at ε =");
  auto unordered         = semicircle(D);
  unordered.energy[10]   = unordered.energy[9];
  expect_invalid_argument([&] { tabulated_lattice_t{unordered}; }, "not strictly increasing at point 10");

  auto bad                  = make_sigma(mesh, sigma_up);
  bad[1].data()(3, 0, 0) = {std::numeric_limits<double>::quiet_NaN(), 0.0};
  expect_invalid_argument([&] { bethe_lattice_t(D).moments(bad, mu); }, "Sigma is not finite in block up, element (0, 0)");

  // A defective Σ (a Jordan block) has no well-conditioned eigenvectors
  auto jordan = make_sigma(mesh, [](double) { return make_matrix({{dcomplex(0, -0.1), 1.0}, {0.0, dcomplex(0, -0.1)}}); });
  expect_invalid_argument([&] { bethe_lattice_t(D).moments(jordan, mu); }, "in block up at ω =");
}
