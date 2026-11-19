#include "./causality.hpp"

#include "../nrg/kramers_kronig.hpp"
#include "../nrg/convert.hpp"

#include <nda/linalg.hpp>

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace nrgljubljana_interface::dmft {

  namespace {

    using dcomplex = std::complex<double>;
    using values_t = nda::array<dcomplex, 3>; // (mesh point, i, j)

    // An eigenvalue counts as below the floor only if it is below it by more than this times the largest |F_ij(ω)|
    constexpr double floor_tolerance = 1e-14;

    // A high-frequency limit is Hermitian if |L_ij − conj(L_ji)| is at most this times the largest |L_ij|
    constexpr double hermitian_tolerance = 1e-12;

    std::size_t workers() { return static_cast<std::size_t>(std::max(1, omp_get_max_threads())); }

    using detail::check_finite;
    using detail::number;

    // The Hermitian part (F + F†)/2 and the spectral part (F† − F)/(2i) of one block, both exactly Hermitian
    void split(values_t const &d, values_t &H, values_t &S) {
      long const nw = d.extent(0), n = d.extent(1);
      H = values_t(nw, n, n);
      S = values_t(nw, n, n);
      for (long k = 0; k < nw; ++k)
        for (long i = 0; i < n; ++i)
          for (long j = 0; j < n; ++j) {
            dcomplex const sum  = d(k, i, j) + std::conj(d(k, j, i));
            dcomplex const diff = std::conj(d(k, j, i)) - d(k, i, j);
            H(k, i, j)          = sum / 2.0;
            S(k, i, j)          = {diff.imag() / 2.0, -diff.real() / 2.0}; // diff / (2i)
          }
    }

    // F = H − iS, element by element
    dcomplex combine(dcomplex h, dcomplex s) { return {h.real() + s.imag(), h.imag() - s.real()}; }

    // KK(−X) of a Hermitian matrix function on the mesh points: the upper triangle element by element (real and
    // imaginary parts separately), the lower triangle by conjugation. Parts that are zero everywhere give zero.
    values_t kk_minus(values_t const &X, std::vector<double> const &points) {
      long const nw = X.extent(0), n = X.extent(1);
      values_t K(nw, n, n);
      K = 0;
      std::vector<double> re(nw), im(nw);
      for (long i = 0; i < n; ++i)
        for (long j = i; j < n; ++j) {
          bool any_re = false, any_im = false;
          for (long k = 0; k < nw; ++k) {
            re[k] = -X(k, i, j).real();
            im[k] = -X(k, i, j).imag();
            any_re = any_re or re[k] != 0;
            any_im = any_im or (i != j and im[k] != 0);
          }
          std::vector<double> const kre = any_re ? nrg::kramers_kronig(points, re, workers()) : std::vector<double>(nw, 0.0);
          std::vector<double> const kim = any_im ? nrg::kramers_kronig(points, im, workers()) : std::vector<double>(nw, 0.0);
          for (long k = 0; k < nw; ++k) {
            K(k, i, j) = {kre[k], kim[k]};
            if (i != j) K(k, j, i) = std::conj(K(k, i, j));
          }
        }
      return K;
    }

    double largest_abs(values_t const &d, long k) {
      double m = 0;
      for (long i = 0; i < d.extent(1); ++i)
        for (long j = 0; j < d.extent(2); ++j) m = std::max(m, std::abs(d(k, i, j)));
      return m;
    }

    // The high-frequency limit checked against F and symmetrized
    block_matrices_t checked_limit(g_w_t const &F, block_matrices_t const &limit) {
      detail::check_block_matrices(F, limit, "rebuild_hermitian_part", "the high-frequency limit");
      block_matrices_t result;
      for (long b = 0; b < long(F.size()); ++b) {
        auto const &L = limit[b];
        long const n  = L.extent(0);
        double scale = 0, asymmetry = 0;
        for (long i = 0; i < n; ++i)
          for (long j = 0; j < n; ++j) {
            scale     = std::max(scale, std::abs(L(i, j)));
            asymmetry = std::max(asymmetry, std::abs(L(i, j) - std::conj(L(j, i))));
          }
        if (asymmetry > hermitian_tolerance * scale)
          throw std::invalid_argument("rebuild_hermitian_part: the high-frequency limit of block " + F.block_names()[b]
                                      + " is not Hermitian (|L - L†| = " + number(asymmetry) + ")");
        nda::matrix<dcomplex> sym(n, n);
        for (long i = 0; i < n; ++i)
          for (long j = 0; j < n; ++j) sym(i, j) = (L(i, j) + std::conj(L(j, i))) / 2.0;
        result.push_back(std::move(sym));
      }
      return result;
    }

  } // namespace

  g_w_t spectral_part(g_w_t const &F) {
    g_w_t S = F;
    for (long b = 0; b < long(F.size()); ++b) {
      values_t H, s;
      split(F[b].data(), H, s);
      S[b].data() = s;
    }
    return S;
  }

  std::pair<g_w_t, causality_report_t> causal_projection(g_w_t const &F, double floor) {
    if (not std::isfinite(floor) or floor < 0) throw std::invalid_argument("causal_projection: the floor must be finite and >= 0, got " + number(floor));
    check_finite(F, "causal_projection");

    g_w_t result = F;
    causality_report_t report;
    for (long b = 0; b < long(F.size()); ++b) {
      auto const &d = F[b].data();
      long const nw = d.extent(0), n = d.extent(1);
      values_t H, S;
      split(d, H, S);
      values_t Sp = S, dS(nw, n, n);
      dS          = 0;
      bool raised = false;

      for (long k = 0; k < nw; ++k) {
        double const limit = floor - floor_tolerance * largest_abs(d, k);
        if (n == 1) {
          double const lambda = S(k, 0, 0).real();
          if (not(lambda < limit)) continue;
          Sp(k, 0, 0) = floor;
          dS(k, 0, 0) = floor - lambda;
          raised      = true;
          report.projected_points++;
          report.largest_violation = std::max(report.largest_violation, floor - lambda);
          continue;
        }
        nda::matrix<dcomplex> s(n, n);
        for (long i = 0; i < n; ++i)
          for (long j = 0; j < n; ++j) s(i, j) = S(k, i, j);
        auto const [lambda, V] = nda::linalg::eigh(s);
        if (not(lambda(0) < limit)) continue;
        // δS = Σ_m (floor − λ_m) v_m v_m† over the raised eigenvalues; upper triangle, the rest by conjugation
        for (long i = 0; i < n; ++i)
          for (long j = i; j < n; ++j) {
            dcomplex sum = 0;
            for (long m = 0; m < n; ++m)
              if (lambda(m) < limit) sum += (floor - lambda(m)) * V(i, m) * std::conj(V(j, m));
            if (i == j) sum = {sum.real(), 0.0};
            dS(k, i, j) = sum;
            if (i != j) dS(k, j, i) = std::conj(sum);
          }
        for (long i = 0; i < n; ++i)
          for (long j = 0; j < n; ++j) Sp(k, i, j) = S(k, i, j) + dS(k, i, j);
        raised = true;
        report.projected_points++;
        report.largest_violation = std::max(report.largest_violation, floor - lambda(0));
      }
      if (not raised) continue;

      auto const K = kk_minus(dS, nrg::mesh_points(F[b].mesh()));
      auto &out    = result[b].data();
      for (long k = 0; k < nw; ++k)
        for (long i = 0; i < n; ++i)
          for (long j = 0; j < n; ++j) out(k, i, j) = combine(H(k, i, j) + K(k, i, j), Sp(k, i, j));
    }
    return {std::move(result), report};
  }

  g_w_t rebuild_hermitian_part(g_w_t const &F, block_matrices_t const &high_frequency_limit) {
    check_finite(F, "rebuild_hermitian_part");
    auto const limit = checked_limit(F, high_frequency_limit);

    g_w_t result = F;
    for (long b = 0; b < long(F.size()); ++b) {
      auto const &d = F[b].data();
      long const nw = d.extent(0), n = d.extent(1);
      values_t H, S;
      split(d, H, S);
      auto const K = kk_minus(S, nrg::mesh_points(F[b].mesh()));
      auto &out    = result[b].data();
      for (long k = 0; k < nw; ++k)
        for (long i = 0; i < n; ++i)
          for (long j = 0; j < n; ++j) out(k, i, j) = combine(limit[b](i, j) + K(k, i, j), S(k, i, j));
    }
    return result;
  }

  block_matrices_t high_frequency_limit(g_w_t const &F) {
    check_finite(F, "high_frequency_limit");
    block_matrices_t result;
    for (long b = 0; b < long(F.size()); ++b) {
      auto const &d = F[b].data();
      long const nw = d.extent(0), n = d.extent(1);
      values_t H, S;
      split(d, H, S);
      auto const K = kk_minus(S, nrg::mesh_points(F[b].mesh()));
      nda::matrix<dcomplex> L(n, n);
      for (long i = 0; i < n; ++i)
        for (long j = 0; j < n; ++j) L(i, j) = ((H(0, i, j) - K(0, i, j)) + (H(nw - 1, i, j) - K(nw - 1, i, j))) / 2.0;
      nda::matrix<dcomplex> sym(n, n);
      for (long i = 0; i < n; ++i)
        for (long j = 0; j < n; ++j) sym(i, j) = (L(i, j) + std::conj(L(j, i))) / 2.0;
      result.push_back(std::move(sym));
    }
    return result;
  }

} // namespace nrgljubljana_interface::dmft
