// Causality kernels (dmft/kernels/causality). No NRG run.
//
// The analytic test functions are bands with compact support, B(ω) = H_B(ω) − i S_B(ω) with
// S_B(ω) = G0 (1 − ((ω − c)/D)²)² for |ω − c| < D, and the closed form of H_B = KK(−S_B). Inside the mesh there is
// no tail, so the Kramers-Kronig transform of the sampled S_B is limited only by the interpolation. The tolerances
// below were measured with the kk tool of NRG Ljubljana (-a analytic -i steffen, the algorithm of the library) on
// the same functions and meshes.

#include <gtest/gtest.h>

#include <nrgljubljana_interface/dmft/kernels/causality.hpp>
#include <nrgljubljana_interface/nrg/kramers_kronig.hpp>
#include <nrgljubljana_interface/nrg/convert.hpp>

#include <nda/linalg.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace nrgljubljana_interface;
using namespace nrgljubljana_interface::dmft;
using triqs::mesh::refreq_log;
using dcomplex = std::complex<double>;
using cmatrix  = nda::matrix<dcomplex>;

namespace {

  // Few points, for the cases that check structure rather than accuracy
  refreq_log small_mesh() { return {1e-2, 5.0, 1.2}; }

  // The accuracy mesh: 780 points
  refreq_log fine_mesh() { return {1e-4, 10.0, 1.03}; }

  struct band_t {
    double c, D, G0;
  };

  dcomplex band(double w, band_t const &b) {
    double const u   = w - b.c;
    double const P   = b.G0 * std::pow(1 - (u / b.D) * (u / b.D), 2);
    double const lg  = std::abs(std::abs(u) - b.D) == 0 ? 0.0 : std::log(std::abs((u + b.D) / (u - b.D)));
    double const H   = (P * lg + b.G0 * (10 * u / (3 * b.D) - 2 * u * u * u / (b.D * b.D * b.D))) / std::numbers::pi;
    double const S_b = std::abs(u) < b.D ? P : 0.0;
    return {H, -S_b};
  }

  // F(ω) = F_∞ + Σ_m M_m B_m(ω), one block
  struct term_t {
    cmatrix M;
    band_t b;
  };

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

  // The matrix at mesh point k of a block's data, and back
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

  g_w_t make_gf(refreq_log const &mesh, std::vector<std::string> const &names, std::vector<cmatrix> const &limits,
                std::vector<std::vector<term_t>> const &terms) {
    triqs::gfs::gf_struct_t s;
    for (std::size_t b = 0; b < names.size(); ++b) s.emplace_back(names[b], limits[b].extent(0));
    g_w_t F{mesh, s};
    for (std::size_t b = 0; b < names.size(); ++b) {
      auto &d         = F[long(b)].data();
      auto const &pts = mesh.points();
      for (long k = 0; k < long(pts.size()); ++k) {
        cmatrix v = limits[b];
        for (auto const &t : terms[b]) v += t.M * band(pts[k], t.b);
        put(d, k, v);
      }
    }
    return F;
  }

  // The three cases with violations
  cmatrix one(double x) { return make_matrix({{x}}); }

  g_w_t plateau(refreq_log const &mesh) { return make_gf(mesh, {"imp"}, {one(0.25)}, {{{one(1.0), {0.3, 1.0, 0.7}}}}); }

  g_w_t notch(refreq_log const &mesh) {
    return make_gf(mesh, {"imp"}, {one(0.25)}, {{{one(1.0), {0.0, 1.0, 0.7}}, {one(-1.0), {0.5, 0.15, 0.42}}}});
  }

  cmatrix limit_2x2() { return make_matrix({{0.1, {0.05, -0.02}}, {{0.05, 0.02}, -0.2}}); }
  cmatrix M1() { return make_matrix({{1.0, 0.0}, {0.0, 0.3}}); }
  cmatrix M2() { return make_matrix({{0.3, {0.0, 0.2}}, {{0.0, -0.2}, 1.0}}); }
  cmatrix M3() {
    dcomplex const v0 = 1.0, v1 = {0.5, 0.5};
    double const norm = std::sqrt(std::norm(v0) + std::norm(v1));
    return make_matrix({{-std::norm(v0) / (norm * norm), -v0 * std::conj(v1) / (norm * norm)},
                   {-v1 * std::conj(v0) / (norm * norm), -std::norm(v1) / (norm * norm)}});
  }

  // Positive semidefinite spectral part everywhere
  std::vector<term_t> causal_2x2() { return {{M1(), {-0.3, 1.0, 0.6}}, {M2(), {0.4, 0.8, 0.5}}}; }

  // A negative eigenvalue near ω = 0.1, with eigenvectors that rotate with ω
  std::vector<term_t> rotating_2x2() {
    auto t = causal_2x2();
    t.push_back({M3(), {0.1, 0.2, 0.5}});
    return t;
  }

  g_w_t rotating(refreq_log const &mesh) { return make_gf(mesh, {"up"}, {limit_2x2()}, {rotating_2x2()}); }

  double max_abs_diff(g_w_t const &a, g_w_t const &b) {
    double m = 0;
    for (long bl = 0; bl < long(a.size()); ++bl) m = std::max(m, double(max_element(abs(a[bl].data() - b[bl].data()))));
    return m;
  }

  bool bitwise_equal(g_w_t const &a, g_w_t const &b) {
    for (long bl = 0; bl < long(a.size()); ++bl) {
      auto const &x = a[bl].data();
      auto const &y = b[bl].data();
      for (long k = 0; k < x.extent(0); ++k)
        for (long i = 0; i < x.extent(1); ++i)
          for (long j = 0; j < x.extent(2); ++j)
            if (x(k, i, j).real() != y(k, i, j).real() or x(k, i, j).imag() != y(k, i, j).imag()) return false;
    }
    return true;
  }

  // The Hermitian part (F + F†)/2
  g_w_t hermitian_part(g_w_t const &F) {
    g_w_t H = F;
    for (long bl = 0; bl < long(F.size()); ++bl) {
      auto const &d = F[bl].data();
      auto &h       = H[bl].data();
      for (long k = 0; k < d.extent(0); ++k)
        for (long i = 0; i < d.extent(1); ++i)
          for (long j = 0; j < d.extent(2); ++j) h(k, i, j) = (d(k, i, j) + std::conj(d(k, j, i))) / 2.0;
    }
    return H;
  }

  // KK(−X) of one Hermitian block, element by element
  nda::array<dcomplex, 3> kk_minus(nda::array<dcomplex, 3> const &X, std::vector<double> const &points) {
    long const nw = X.extent(0), n = X.extent(1);
    nda::array<dcomplex, 3> K(nw, n, n);
    std::vector<double> re(nw), im(nw);
    for (long i = 0; i < n; ++i)
      for (long j = 0; j < n; ++j) {
        for (long k = 0; k < nw; ++k) {
          re[k] = -X(k, i, j).real();
          im[k] = -X(k, i, j).imag();
        }
        auto const kre = nrg::kramers_kronig(points, re);
        auto const kim = nrg::kramers_kronig(points, im);
        for (long k = 0; k < nw; ++k) K(k, i, j) = {kre[k], kim[k]};
      }
    return K;
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

TEST(Causality, SpectralPart) {
  auto const mesh = small_mesh();
  g_w_t F{mesh, triqs::gfs::gf_struct_t{{"a", 1}, {"b", 2}}};
  auto const &pts = mesh.points();
  for (long k = 0; k < long(pts.size()); ++k) {
    F[0].data()(k, 0, 0) = {0.3 + k, -0.7 - 0.01 * k};
    for (long i = 0; i < 2; ++i)
      for (long j = 0; j < 2; ++j) F[1].data()(k, i, j) = {1.0 + 0.1 * i + 0.01 * j + k, -2.0 + 0.3 * i - 0.07 * j - 0.5 * k};
  }
  auto const S = spectral_part(F);
  for (long k = 0; k < long(pts.size()); ++k) {
    EXPECT_EQ(S[0].data()(k, 0, 0), dcomplex(-F[0].data()(k, 0, 0).imag(), 0.0));
    for (long i = 0; i < 2; ++i)
      for (long j = 0; j < 2; ++j) {
        dcomplex const expected = (std::conj(F[1].data()(k, j, i)) - F[1].data()(k, i, j)) / dcomplex(0, 2);
        EXPECT_LT(std::abs(S[1].data()(k, i, j) - expected), 1e-14 * (1 + k));
        EXPECT_EQ(S[1].data()(k, j, i), std::conj(S[1].data()(k, i, j)));
      }
  }
}

TEST(Causality, NoViolationUnchanged) {
  auto const mesh = small_mesh();
  auto const F    = make_gf(mesh, {"imp", "up"}, {one(0.25), limit_2x2()}, {{{one(1.0), {0.3, 1.0, 0.7}}}, causal_2x2()});
  auto const [P, report] = causal_projection(F, 0.0);
  EXPECT_TRUE(bitwise_equal(P, F));
  EXPECT_EQ(report, causality_report_t{});
}

TEST(Causality, ScalarViolation) {
  auto const mesh   = small_mesh();
  auto const F      = plateau(mesh);
  double const fl   = 1e-3;
  auto const &d     = F[0].data();
  auto const points = nrg::mesh_points(mesh);
  long const nw     = d.extent(0);

  std::vector<double> dS(nw, 0.0);
  causality_report_t expected;
  for (long k = 0; k < nw; ++k) {
    double const s = -d(k, 0, 0).imag();
    if (s < fl - 1e-14 * std::abs(d(k, 0, 0))) {
      dS[k] = fl - s;
      expected.projected_points++;
      expected.largest_violation = std::max(expected.largest_violation, fl - s);
    }
  }
  ASSERT_GT(expected.projected_points, 0);
  std::vector<double> minus(nw);
  for (long k = 0; k < nw; ++k) minus[k] = -dS[k];
  auto const K = nrg::kramers_kronig(points, minus);

  auto const [P, report] = causal_projection(F, fl);
  EXPECT_EQ(report, expected);
  auto const &p = P[0].data();
  for (long k = 0; k < nw; ++k) {
    EXPECT_EQ(p(k, 0, 0).imag(), dS[k] > 0 ? -fl : d(k, 0, 0).imag()) << "k = " << k;
    EXPECT_EQ(p(k, 0, 0).real(), d(k, 0, 0).real() + K[k]) << "k = " << k;
  }
}

TEST(Causality, MatrixViolation) {
  auto const mesh        = small_mesh();
  auto const F           = rotating(mesh);
  auto const [P, report] = causal_projection(F, 0.0);
  ASSERT_GT(report.projected_points, 0);
  EXPECT_GT(report.largest_violation, 1e-3);

  auto const S  = spectral_part(F)[0].data();
  auto const Sp = spectral_part(P)[0].data();
  long const nw = S.extent(0);
  double scale  = 0;
  for (long k = 0; k < nw; ++k) scale = std::max(scale, double(max_element(abs(at(F[0].data(), k)))));
  for (long k = 0; k < nw; ++k) {
    cmatrix const s  = at(S, k);
    cmatrix const sp = at(Sp, k);
    auto const lambda = nda::linalg::eigvalsh(sp);
    EXPECT_GE(lambda(0), -1e-13 * scale) << "k = " << k;
    cmatrix const commutator = sp * s - s * sp;
    EXPECT_LT(double(max_element(abs(commutator))), 1e-12 * scale * scale) << "k = " << k;
  }
  // The Hermitian part changed by KK(−δS)
  nda::array<dcomplex, 3> const dS = Sp - S;
  auto const K               = kk_minus(dS, nrg::mesh_points(mesh));
  nda::array<dcomplex, 3> const change = hermitian_part(P)[0].data() - hermitian_part(F)[0].data();
  EXPECT_LT(double(max_element(abs(change - K))), 1e-12 * scale);
}

// For an analytic F, the projection equals the rebuild of the projected F with the true F_∞, up to the Kramers-Kronig
// error the case has anyway (measured by rebuilding the unprojected F). Measured ratios of the two errors: plateau
// 1.0, notch 1.3, rotating 2x2 2.8. The correction stays of the order of the change of S.
TEST(Causality, ProjectionMatchesRebuild) {
  auto const mesh = fine_mesh();
  struct case_t {
    std::string name;
    g_w_t F;
    cmatrix limit;
    double floor;
    long projected_points;
  };
  std::vector<case_t> cases{{"plateau", plateau(mesh), one(0.25), 1e-6, 160},
                            {"notch", notch(mesh), one(0.25), 0.0, 4},
                            {"rotating 2x2", rotating(mesh), limit_2x2(), 0.0, 26}};
  for (auto const &c : cases) {
    auto const base        = max_abs_diff(rebuild_hermitian_part(c.F, {c.limit}), c.F);
    auto const [P, report] = causal_projection(c.F, c.floor);
    EXPECT_EQ(report.projected_points, c.projected_points) << c.name;
    auto const rebuilt = rebuild_hermitian_part(P, {c.limit});
    EXPECT_LT(max_abs_diff(P, rebuilt), 4 * base) << c.name << ", rebuild error of the unprojected function " << base;

    nda::array<dcomplex, 3> const dS     = spectral_part(P)[0].data() - spectral_part(c.F)[0].data();
    nda::array<dcomplex, 3> const change = hermitian_part(P)[0].data() - hermitian_part(c.F)[0].data();
    EXPECT_LT(double(max_element(abs(change))), 5 * double(max_element(abs(dS)))) << c.name;
  }
}

// Measured on the fine mesh: 3.3e-5 for the 1x1 band, 2.7e-4 for the 2x2 function (at the edges of the 0.8-wide band)
TEST(Causality, RebuildAnalytic) {
  auto const mesh = fine_mesh();
  auto const F1   = make_gf(mesh, {"imp"}, {one(0.25)}, {{{one(1.0), {0.3, 1.0, 0.7}}}});
  EXPECT_LT(max_abs_diff(rebuild_hermitian_part(F1, {one(0.25)}), F1), 1e-4);
  auto const F2 = make_gf(mesh, {"up"}, {limit_2x2()}, {causal_2x2()});
  EXPECT_LT(max_abs_diff(rebuild_hermitian_part(F2, {limit_2x2()}), F2), 5e-4);
}

// Measured on the fine mesh: errors 3e-9 (1x1) and 7e-9 (2x2); the plain average of H over the two ends is off by
// 7e-4 and 4.5e-4, the first moment of ρ over ω_max²
TEST(Causality, HighFrequencyLimit) {
  auto const mesh = fine_mesh();
  auto const F    = make_gf(mesh, {"imp", "up"}, {one(0.25), limit_2x2()}, {{{one(1.0), {0.3, 1.0, 0.7}}}, causal_2x2()});
  auto const L    = high_frequency_limit(F);
  ASSERT_EQ(L.size(), 2u);
  std::vector<cmatrix> const exact{one(0.25), limit_2x2()};
  auto const H = hermitian_part(F);
  for (long b = 0; b < 2; ++b) {
    auto const &h     = H[b].data();
    long const nw     = h.extent(0);
    cmatrix const plain = (at(h, 0) + at(h, nw - 1)) / 2.0;
    double const err  = double(max_element(abs(L[b] - exact[b])));
    double const perr = double(max_element(abs(plain - exact[b])));
    EXPECT_LT(err, 1e-7) << "block " << b;
    EXPECT_GT(perr, 1e-4) << "block " << b;
    if (b == 1) { EXPECT_EQ(L[b](0, 1), std::conj(L[b](1, 0))); }
  }
}

// A unitary rotation of the channels: the spectral parts and the report rotate exactly (up to rounding). The Hermitian
// parts come from KK, which interpolates each element with Steffen; that is not linear in the data, so they rotate
// only to the Kramers-Kronig error the case has anyway (the rebuild error of the unprojected function).
TEST(Causality, UnitaryRotation) {
  auto const mesh  = fine_mesh();
  auto const F     = rotating(mesh);
  cmatrix const L  = limit_2x2();
  double const t   = 0.7, p = 0.4;
  cmatrix const U  = make_matrix({{std::cos(t), -std::sin(t) * std::exp(dcomplex(0, -p))}, {std::sin(t) * std::exp(dcomplex(0, p)), std::cos(t)}});
  cmatrix const Ud = nda::dagger(U);
  auto rotate      = [&](g_w_t const &G) {
    g_w_t R = G;
    for (long k = 0; k < G[0].data().extent(0); ++k) put(R[0].data(), k, cmatrix(U * at(G[0].data(), k) * Ud));
    return R;
  };
  auto const base = max_abs_diff(rebuild_hermitian_part(F, {L}), F);

  auto const [P, report]    = causal_projection(F, 0.0);
  auto const [PR, report_r] = causal_projection(rotate(F), 0.0);
  EXPECT_EQ(report.projected_points, report_r.projected_points);
  EXPECT_NEAR(report.largest_violation, report_r.largest_violation, 1e-13);
  EXPECT_LT(max_abs_diff(spectral_part(PR), rotate(spectral_part(P))), 1e-12);
  EXPECT_LT(max_abs_diff(PR, rotate(P)), base) << "rebuild error of the unprojected function " << base;

  auto const B  = rebuild_hermitian_part(F, {L});
  auto const BR = rebuild_hermitian_part(rotate(F), {cmatrix(U * L * Ud)});
  EXPECT_LT(max_abs_diff(spectral_part(BR), rotate(spectral_part(B))), 1e-12);
  EXPECT_LT(max_abs_diff(BR, rotate(B)), 2 * base) << "rebuild error of the unprojected function " << base;
}

TEST(Causality, SecondProjectionChangesNothing) {
  auto const mesh = small_mesh();
  for (auto const &[F, fl] : std::vector<std::pair<g_w_t, double>>{{plateau(mesh), 1e-3}, {rotating(mesh), 0.0}, {rotating(mesh), 1e-3}}) {
    auto const first           = causal_projection(F, fl).first;
    auto const [second, again] = causal_projection(first, fl);
    EXPECT_EQ(again, causality_report_t{});
    EXPECT_TRUE(bitwise_equal(second, first));
  }
}

TEST(Causality, Errors) {
  auto const mesh = small_mesh();
  auto const F    = make_gf(mesh, {"imp", "up"}, {one(0.25), limit_2x2()}, {{{one(1.0), {0.3, 1.0, 0.7}}}, causal_2x2()});

  expect_invalid_argument([&] { causal_projection(F, -1.0); }, "floor must be finite and >= 0");
  expect_invalid_argument([&] { causal_projection(F, std::numeric_limits<double>::quiet_NaN()); }, "floor must be finite and >= 0");

  auto bad                  = F;
  bad[1].data()(3, 0, 1) = {std::numeric_limits<double>::infinity(), 0.0};
  expect_invalid_argument([&] { causal_projection(bad, 0.0); }, "not finite in block up, element (0, 1)");
  expect_invalid_argument([&] { rebuild_hermitian_part(bad, {one(0.25), limit_2x2()}); }, "not finite in block up");
  expect_invalid_argument([&] { high_frequency_limit(bad); }, "not finite in block up");

  expect_invalid_argument([&] { rebuild_hermitian_part(F, {one(0.25)}); }, "1 matrices, the function 2 blocks");
  expect_invalid_argument([&] { rebuild_hermitian_part(F, {one(0.25), one(0.0)}); }, "block up is 1x1, the block 2x2");
  expect_invalid_argument([&] { rebuild_hermitian_part(F, {one(0.25), make_matrix({{0.1, 0.2}, {0.3, 0.0}})}); }, "block up is not Hermitian");
  expect_invalid_argument([&] { rebuild_hermitian_part(F, {make_matrix({{dcomplex(0.25, 1e-3)}}), limit_2x2()}); }, "block imp is not Hermitian");
  expect_invalid_argument([&] { rebuild_hermitian_part(F, {one(std::numeric_limits<double>::quiet_NaN()), limit_2x2()}); },
                          "block imp is not finite");
}

TEST(Causality, KramersKronigWorkers) {
  auto const mesh   = fine_mesh();
  auto const points = nrg::mesh_points(mesh);
  std::vector<double> y(points.size());
  for (std::size_t k = 0; k < points.size(); ++k) y[k] = band(points[k], {0.3, 1.0, 0.7}).imag();
  EXPECT_EQ(nrg::kramers_kronig(points, y, 1), nrg::kramers_kronig(points, y, 4));
}
