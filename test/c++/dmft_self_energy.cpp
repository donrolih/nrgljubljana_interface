// Self-energy functions (dmft/self_energy). No NRG run.
//
// Constructed case: a self-energy Σ = Σ_H + Σ', a hybridization Δ and levels h, G = [ω − h − Δ − Σ]⁻¹, and the
// correlators F_l = Σ'G, F_r = GΣ', I = Σ' + Σ'GΣ', for which Σ_H + I − F_l G⁻¹ F_r = Σ and F_l G⁻¹ = Σ' exactly.
// Blocks of size 1 and 2, with matrices that do not commute.

#include <gtest/gtest.h>

#include <nrgljubljana_interface/dmft/self_energy.hpp>

#include <nda/linalg.hpp>

#include <algorithm>
#include <complex>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace nrgljubljana_interface;
using namespace nrgljubljana_interface::dmft;
using triqs::mesh::refreq_log;
using dcomplex = std::complex<double>;
using cmatrix  = nda::matrix<dcomplex>;

namespace {

  refreq_log test_mesh() { return {1e-2, 5.0, 1.2}; }

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

  void put(nda::array<dcomplex, 3> &d, long k, cmatrix const &m) {
    for (long i = 0; i < d.extent(1); ++i)
      for (long j = 0; j < d.extent(2); ++j) d(k, i, j) = m(i, j);
  }

  // A Green's function with blocks imp (1x1) and up (2x2), given per block as a function of ω
  template <typename Function> g_w_t make_gf(refreq_log const &mesh, Function f) {
    g_w_t g{mesh, triqs::gfs::gf_struct_t{{"imp", 1}, {"up", 2}}};
    auto const &pts = mesh.points();
    for (long b = 0; b < 2; ++b)
      for (long k = 0; k < long(pts.size()); ++k) put(g[b].data(), k, f(b, pts[k]));
    return g;
  }

  // A Lorentzian w/(ω − c + i w), causal
  dcomplex lorentz(double w, double c, double width) { return width / (w - c + dcomplex(0, width)); }

  cmatrix id(long b) { return b == 0 ? make_matrix({{1.0}}) : make_matrix({{1.0, 0.0}, {0.0, 1.0}}); }

  // Σ', Δ and the constants of the constructed case
  cmatrix sigma_prime(long b, double w) {
    if (b == 0) return make_matrix({{0.8 * lorentz(w, 0.3, 0.5) - dcomplex(0, 0.1)}});
    return make_matrix({{1.0, {0.2, 0.3}}, {{0.2, -0.3}, 0.5}}) * lorentz(w, -0.4, 0.6)
       + make_matrix({{0.3, 0.1}, {0.1, 0.9}}) * lorentz(w, 0.5, 0.3) - dcomplex(0, 0.1) * id(1);
  }

  cmatrix delta(long b, double w) {
    if (b == 0) return make_matrix({{0.5 * lorentz(w, 0.0, 1.0)}});
    return make_matrix({{0.6, {0.0, 0.2}}, {{0.0, -0.2}, 0.4}}) * lorentz(w, 0.1, 0.8);
  }

  block_matrices_t hartree() { return {make_matrix({{1.7}}), make_matrix({{1.1, 0.2}, {0.2, 0.9}})}; }
  block_matrices_t levels() { return {make_matrix({{-0.5}}), make_matrix({{-0.3, {0.05, 0.1}}, {{0.05, -0.1}, -0.6}})}; }

  cmatrix g_exact(long b, double w) { return nda::linalg::inv(cmatrix(w * id(b) - levels()[b] - delta(b, w) - hartree()[b] - sigma_prime(b, w))); }

  struct constructed_t {
    g_w_t G, F_l, F_r, I, Delta, Sigma, Sigma_prime;
  };

  constructed_t constructed(refreq_log const &mesh) {
    auto G   = make_gf(mesh, g_exact);
    auto F_l = make_gf(mesh, [](long b, double w) { return cmatrix(sigma_prime(b, w) * g_exact(b, w)); });
    auto F_r = make_gf(mesh, [](long b, double w) { return cmatrix(g_exact(b, w) * sigma_prime(b, w)); });
    auto I   = make_gf(mesh, [](long b, double w) { return cmatrix(sigma_prime(b, w) + sigma_prime(b, w) * g_exact(b, w) * sigma_prime(b, w)); });
    auto D   = make_gf(mesh, delta);
    auto S   = make_gf(mesh, [](long b, double w) { return cmatrix(hartree()[b] + sigma_prime(b, w)); });
    auto Sp  = make_gf(mesh, sigma_prime);
    return {G, F_l, F_r, I, D, S, Sp};
  }

  // Largest |a − b| relative to the largest |b|
  double relative_difference(g_w_t const &a, g_w_t const &b) {
    double diff = 0, scale = 0;
    for (long bl = 0; bl < long(a.size()); ++bl) {
      diff  = std::max(diff, double(max_element(abs(a[bl].data() - b[bl].data()))));
      scale = std::max(scale, double(max_element(abs(b[bl].data()))));
    }
    return diff / scale;
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

} // namespace

TEST(SelfEnergy, SigmaIFG) {
  auto const c = constructed(test_mesh());
  EXPECT_LT(relative_difference(sigma_IFG(c.G, c.F_l, c.F_r, c.I, hartree()), c.Sigma), 1e-13);
}

TEST(SelfEnergy, SigmaFG) {
  auto const c = constructed(test_mesh());
  EXPECT_LT(relative_difference(sigma_FG(c.G, c.F_l), c.Sigma_prime), 1e-13);
}

TEST(SelfEnergy, ImpurityGf) {
  auto const c = constructed(test_mesh());
  auto const G = impurity_gf(c.Delta, c.Sigma, levels());
  EXPECT_LT(relative_difference(G, c.G), 1e-13);
  // Its inverse is ω − h − Δ − Σ
  auto const mesh = test_mesh();
  auto const &pts = mesh.points();
  for (long b = 0; b < 2; ++b)
    for (long k = 0; k < long(pts.size()); ++k) {
      cmatrix const expected = pts[k] * id(b) - levels()[b] - at(c.Delta[b].data(), k) - at(c.Sigma[b].data(), k);
      cmatrix const inverse  = nda::linalg::inv(at(G[b].data(), k));
      EXPECT_LT(double(max_element(abs(inverse - expected))), 1e-12 * (1 + double(max_element(abs(expected)))));
    }
}

TEST(SelfEnergy, Errors) {
  auto const mesh = test_mesh();
  auto const c    = constructed(mesh);

  g_w_t other_blocks{mesh, triqs::gfs::gf_struct_t{{"imp", 1}, {"dn", 2}}};
  expect_invalid_argument([&] { sigma_IFG(c.G, other_blocks, c.F_r, c.I, hartree()); }, "sigma_IFG: F_l does not have the blocks of G");
  g_w_t other_size{mesh, triqs::gfs::gf_struct_t{{"imp", 1}, {"up", 3}}};
  expect_invalid_argument([&] { sigma_IFG(c.G, c.F_l, other_size, c.I, hartree()); }, "F_r block up is 3x3, in G 2x2");
  g_w_t other_mesh{refreq_log{1e-2, 5.0, 1.3}, triqs::gfs::gf_struct_t{{"imp", 1}, {"up", 2}}};
  expect_invalid_argument([&] { sigma_FG(c.G, other_mesh); }, "F_l block imp is not on the mesh of G");
  expect_invalid_argument([&] { impurity_gf(c.Delta, other_mesh, levels()); }, "Sigma block imp is not on the mesh of Delta");

  expect_invalid_argument([&] { sigma_IFG(c.G, c.F_l, c.F_r, c.I, {hartree()[0]}); }, "the Hartree term has 1 matrices, the function 2 blocks");
  expect_invalid_argument([&] { impurity_gf(c.Delta, c.Sigma, {levels()[0], levels()[0]}); }, "h_imp of block up is 1x1, the block 2x2");
  expect_invalid_argument([&] { impurity_gf(c.Delta, c.Sigma, {make_matrix({{std::numeric_limits<double>::quiet_NaN()}}), levels()[1]}); },
                          "h_imp of block imp is not finite");

  auto bad                  = c.I;
  bad[1].data()(2, 1, 0) = {std::numeric_limits<double>::infinity(), 0.0};
  expect_invalid_argument([&] { sigma_IFG(c.G, c.F_l, c.F_r, bad, hartree()); }, "I is not finite in block up, element (1, 0)");

  auto singular               = c.G;
  singular[0].data()(4, 0, 0) = 0;
  expect_invalid_argument([&] { sigma_IFG(singular, c.F_l, c.F_r, c.I, hartree()); }, "sigma_IFG: G is singular in block imp at ω = ");
  expect_invalid_argument([&] { sigma_FG(singular, c.F_l); }, "sigma_FG: G is singular in block imp");

  // ω − h − Δ − Σ singular at one point: Σ chosen so that it vanishes there
  auto sigma_singular = c.Sigma;
  long const k        = 7;
  double const w      = mesh.points()[k];
  put(sigma_singular[1].data(), k, cmatrix(w * id(1) - levels()[1] - at(c.Delta[1].data(), k)));
  expect_invalid_argument([&] { impurity_gf(c.Delta, sigma_singular, levels()); }, "impurity_gf: ω − h_imp − Δ − Σ is singular in block up");
}
