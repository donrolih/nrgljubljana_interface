#include "./lattice.hpp"

#include <nda/linalg.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace nrgljubljana_interface::dmft {

  namespace {

    using dcomplex = std::complex<double>;
    using cmatrix  = nda::matrix<dcomplex>;

    // Largest condition number of the eigenvector matrix of Σ accepted by the matrix transform
    constexpr double max_condition = 1e8;

    cmatrix at(nda::array<dcomplex, 3> const &d, long k) {
      cmatrix m(d.extent(1), d.extent(2));
      for (long i = 0; i < d.extent(1); ++i)
        for (long j = 0; j < d.extent(2); ++j) m(i, j) = d(k, i, j);
      return m;
    }

    double frobenius(cmatrix const &m) {
      double s = 0;
      for (long i = 0; i < m.extent(0); ++i)
        for (long j = 0; j < m.extent(1); ++j) s += std::norm(m(i, j));
      return std::sqrt(s);
    }

  } // namespace

  g_w_t lattice_t::local_gf(g_w_t const &Sigma, double mu) const { return moments(Sigma, mu).first; }

  std::pair<g_w_t, g_w_t> lattice_t::self_consistency(g_w_t const &Sigma, double mu) const {
    auto [G0, G1] = moments(Sigma, mu);
    g_w_t Delta   = G1;
    for (long b = 0; b < long(G0.size()); ++b) {
      auto const &g0  = G0[b].data();
      auto const &g1  = G1[b].data();
      auto &del       = Delta[b].data();
      auto const &pts = G0[b].mesh().points();
      for (long k = 0; k < g0.extent(0); ++k) {
        cmatrix const m = at(g0, k);
        if (nda::linalg::det(m) == dcomplex{0})
          throw std::invalid_argument("self_consistency: G_loc is singular in block " + G0.block_names()[b] + " at ω = " + detail::number(pts[k]));
        cmatrix const d = at(g1, k) * nda::linalg::inv(m);
        for (long i = 0; i < d.extent(0); ++i)
          for (long j = 0; j < d.extent(1); ++j) del(k, i, j) = d(i, j);
      }
    }
    return {std::move(G0), std::move(Delta)};
  }

  namespace detail {

    void check_Sigma_floor(double floor, std::string const &function) {
      double const smallest = std::sqrt(std::numeric_limits<double>::min());
      if (not std::isfinite(floor) or floor < smallest)
        throw std::invalid_argument(function + ": the floor of Im Σ must be finite and at least " + number(smallest) + ", got " + number(floor));
    }

    std::vector<g_w_t> lattice_moments(g_w_t const &Sigma, double mu, unsigned max_moment, scalar_moment_t const &scalar, bool shared_orbitals,
                                       std::string const &function) {
      if (not std::isfinite(mu)) throw std::invalid_argument(function + ": μ must be finite, got " + number(mu));
      check_finite(Sigma, function, "Sigma");
      std::vector<g_w_t> result(max_moment + 1, Sigma);
      for (auto &G : result)
        for (long b = 0; b < long(G.size()); ++b) G[b].data() = 0;

      for (long b = 0; b < long(Sigma.size()); ++b) {
        auto const &name = Sigma.block_names()[b];
        auto const &d    = Sigma[b].data();
        auto const &mesh = Sigma[b].mesh().points();
        long const nw = d.extent(0), n = d.extent(1);

        double off_diagonal = 0;
        for (long k = 0; k < nw; ++k)
          for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
              if (i != j) off_diagonal = std::max(off_diagonal, std::abs(d(k, i, j)));

        if (off_diagonal == 0) {
          std::vector<dcomplex> sigma(nw);
          for (long i = 0; i < n; ++i) {
            for (long k = 0; k < nw; ++k) sigma[k] = d(k, i, i);
            for (unsigned m = 0; m <= max_moment; ++m) {
              auto const g = scalar(name, i, mesh, sigma, mu, m);
              auto &out    = result[m][b].data();
              for (long k = 0; k < nw; ++k) out(k, i, i) = g[k];
            }
          }
          continue;
        }

        if (not shared_orbitals)
          throw std::invalid_argument(function + ": Σ has off-diagonal elements in block " + name + " (largest " + number(off_diagonal)
                                      + "), which needs one DOS shared by the orbitals of the block");

        // Σ(ω) = U diag(σ_m) U⁻¹ at each point
        std::vector<std::vector<dcomplex>> branches(n, std::vector<dcomplex>(nw));
        std::vector<cmatrix> U(nw), U_inv(nw);
        for (long k = 0; k < nw; ++k) {
          auto [lambda, V] = nda::linalg::eig(at(d, k));
          U[k]             = V;
          if (nda::linalg::det(U[k]) == dcomplex{0})
            throw std::invalid_argument(function + ": Σ is not diagonalizable in block " + name + " at ω = " + number(mesh[k]));
          U_inv[k]               = nda::linalg::inv(U[k]);
          double const condition = frobenius(U[k]) * frobenius(U_inv[k]);
          if (not(condition <= max_condition))
            throw std::invalid_argument(function + ": the eigenvectors of Σ are ill-conditioned in block " + name + " at ω = " + number(mesh[k])
                                        + " (condition number " + number(condition) + ")");
          for (long m = 0; m < n; ++m) branches[m][k] = lambda(m);
        }
        for (unsigned moment = 0; moment <= max_moment; ++moment) {
          std::vector<std::vector<dcomplex>> g(n);
          for (long m = 0; m < n; ++m) g[m] = scalar(name, m, mesh, branches[m], mu, moment);
          auto &out = result[moment][b].data();
          for (long k = 0; k < nw; ++k) {
            cmatrix diag(n, n);
            diag = 0;
            for (long m = 0; m < n; ++m) diag(m, m) = g[m][k];
            cmatrix const G = U[k] * diag * U_inv[k];
            for (long i = 0; i < n; ++i)
              for (long j = 0; j < n; ++j) out(k, i, j) = G(i, j);
          }
        }
      }
      return result;
    }

  } // namespace detail

} // namespace nrgljubljana_interface::dmft
