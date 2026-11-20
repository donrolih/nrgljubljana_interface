#include "./self_energy.hpp"

#include <nda/linalg.hpp>

#include <complex>
#include <stdexcept>
#include <string>

namespace nrgljubljana_interface::dmft {

  namespace {

    using dcomplex = std::complex<double>;

    nda::matrix<dcomplex> at(nda::array<dcomplex, 3> const &d, long k) {
      nda::matrix<dcomplex> m(d.extent(1), d.extent(2));
      for (long i = 0; i < d.extent(1); ++i)
        for (long j = 0; j < d.extent(2); ++j) m(i, j) = d(k, i, j);
      return m;
    }

    // Throws if the matrix of a block of M is singular at some mesh point
    void check_invertible(g_w_t const &M, std::string const &function, std::string const &what) {
      for (long b = 0; b < long(M.size()); ++b) {
        auto const &d   = M[b].data();
        auto const &pts = M[b].mesh().points();
        for (long k = 0; k < d.extent(0); ++k)
          if (nda::linalg::det(at(d, k)) == dcomplex{0})
            throw std::invalid_argument(function + ": " + what + " is singular in block " + M.block_names()[b] + " at ω = " + detail::number(pts[k]));
      }
    }

    g_w_t checked_result(g_w_t result, std::string const &function) {
      detail::check_finite(result, function, "the result");
      return result;
    }

  } // namespace

  g_w_t Sigma_IFG(g_w_t const &G, g_w_t const &F_l, g_w_t const &F_r, g_w_t const &I, block_matrices_t const &Sigma_Hartree) {
    std::string const function = "Sigma_IFG";
    detail::check_same_structure(F_l, G, function, "F_l", "G");
    detail::check_same_structure(F_r, G, function, "F_r", "G");
    detail::check_same_structure(I, G, function, "I", "G");
    detail::check_block_matrices(G, Sigma_Hartree, function, "the Hartree term");
    detail::check_finite(G, function, "G");
    detail::check_finite(F_l, function, "F_l");
    detail::check_finite(F_r, function, "F_r");
    detail::check_finite(I, function, "I");
    check_invertible(G, function, "G");

    g_w_t Sigma = detail::constant_gf(G, Sigma_Hartree) + I - F_l / G * F_r;
    return checked_result(std::move(Sigma), function);
  }

  g_w_t Sigma_FG(g_w_t const &G, g_w_t const &F_l) {
    std::string const function = "Sigma_FG";
    detail::check_same_structure(F_l, G, function, "F_l", "G");
    detail::check_finite(G, function, "G");
    detail::check_finite(F_l, function, "F_l");
    check_invertible(G, function, "G");

    g_w_t Sigma = F_l / G;
    return checked_result(std::move(Sigma), function);
  }

  g_w_t impurity_gf(g_w_t const &Delta, g_w_t const &Sigma, block_matrices_t const &h_imp) {
    std::string const function = "impurity_gf";
    detail::check_same_structure(Sigma, Delta, function, "Sigma", "Delta");
    detail::check_block_matrices(Delta, h_imp, function, "h_imp");
    detail::check_finite(Delta, function, "Delta");
    detail::check_finite(Sigma, function, "Sigma");

    g_w_t inverse = Delta;
    for (long b = 0; b < long(Delta.size()); ++b) {
      auto &d          = inverse[b].data();
      auto const &del  = Delta[b].data();
      auto const &sig  = Sigma[b].data();
      auto const &pts  = Delta[b].mesh().points();
      long const n     = d.extent(1);
      for (long k = 0; k < d.extent(0); ++k)
        for (long i = 0; i < n; ++i)
          for (long j = 0; j < n; ++j) d(k, i, j) = (i == j ? dcomplex{pts[k]} : dcomplex{0}) - h_imp[b](i, j) - del(k, i, j) - sig(k, i, j);
    }
    check_invertible(inverse, function, "ω − h_imp − Δ − Σ");

    g_w_t G = inverse;
    for (long b = 0; b < long(G.size()); ++b) {
      auto &d = G[b].data();
      for (long k = 0; k < d.extent(0); ++k) {
        auto const m = nda::linalg::inv(at(inverse[b].data(), k));
        for (long i = 0; i < d.extent(1); ++i)
          for (long j = 0; j < d.extent(2); ++j) d(k, i, j) = m(i, j);
      }
    }
    return checked_result(std::move(G), function);
  }

} // namespace nrgljubljana_interface::dmft
