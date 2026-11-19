#include "./measurements.hpp"

#include "../nrg/convert.hpp"
#include "../nrg/integrate.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <numbers>
#include <stdexcept>

namespace nrgljubljana_interface::dmft {

  namespace {

    using dcomplex = std::complex<double>;

    // The real and imaginary parts of A_ij = (conj G_ji − G_ij)/(2πi) on the mesh points of one block
    void spectral_element(nda::array<dcomplex, 3> const &d, long i, long j, std::vector<double> &re, std::vector<double> &im) {
      long const nw = d.extent(0);
      re.resize(nw);
      im.resize(nw);
      for (long k = 0; k < nw; ++k) {
        dcomplex const diff = std::conj(d(k, j, i)) - d(k, i, j);
        re[k]               = diff.imag() / (2 * std::numbers::pi);
        im[k]               = -diff.real() / (2 * std::numbers::pi);
      }
    }

    // Per block, a matrix of ∫ of each element of A, given the integral of a real function; Hermitian by construction
    block_matrices_t integrate_elements(g_w_t const &G, std::function<double(std::vector<double> const &, std::vector<double> const &)> const &integral) {
      block_matrices_t result;
      std::vector<double> re, im;
      for (long b = 0; b < long(G.size()); ++b) {
        auto const &d      = G[b].data();
        auto const points  = nrg::mesh_points(G[b].mesh());
        long const n       = d.extent(1);
        nda::matrix<dcomplex> m(n, n);
        for (long i = 0; i < n; ++i)
          for (long j = i; j < n; ++j) {
            spectral_element(d, i, j, re, im);
            double const r = integral(points, re);
            double const c = i == j ? 0.0 : integral(points, im);
            m(i, j)        = {r, c};
            m(j, i)        = {r, -c};
          }
        result.push_back(std::move(m));
      }
      return result;
    }

  } // namespace

  occupancy_t occupancy(g_w_t const &G, double T, std::map<std::string, long> const &degeneracy, double epsabs, double epsrel) {
    if (not std::isfinite(T) or T < 0) throw std::invalid_argument("occupancy: the temperature must be finite and >= 0, got " + detail::number(T));
    detail::check_finite(G, "occupancy", "G");
    auto const &names = G.block_names();
    for (auto const &name : names) {
      auto const it = degeneracy.find(name);
      if (it == degeneracy.end()) throw std::invalid_argument("occupancy: no degeneracy for block " + name);
      if (it->second < 1) throw std::invalid_argument("occupancy: the degeneracy of block " + name + " is " + std::to_string(it->second));
    }
    for (auto const &[name, count] : degeneracy)
      if (std::find(names.begin(), names.end(), name) == names.end()) throw std::invalid_argument("occupancy: degeneracy for block " + name + ", not in G");

    occupancy_t result;
    result.matrix = integrate_elements(G, [&](std::vector<double> const &x, std::vector<double> const &y) {
      return nrg::integrate_fermi(x, y, T, epsabs, epsrel);
    });
    for (long b = 0; b < long(G.size()); ++b) {
      double trace = 0;
      for (long i = 0; i < result.matrix[b].extent(0); ++i) trace += result.matrix[b](i, i).real();
      result.total += double(degeneracy.at(names[b])) * trace;
    }
    return result;
  }

  block_matrices_t spectral_weight(g_w_t const &G) {
    detail::check_finite(G, "spectral_weight", "G");
    return integrate_elements(G, [](std::vector<double> const &x, std::vector<double> const &y) { return nrg::integrate(x, y); });
  }

  distance_t spectral_distance(g_w_t const &a, g_w_t const &b) {
    detail::check_same_structure(b, a, "spectral_distance", "b", "a");
    detail::check_finite(a, "spectral_distance", "a");
    detail::check_finite(b, "spectral_distance", "b");

    distance_t result;
    std::vector<double> re_a, im_a, re_b, im_b, re(0), im(0);
    for (long bl = 0; bl < long(a.size()); ++bl) {
      auto const &da     = a[bl].data();
      auto const &db     = b[bl].data();
      auto const points  = nrg::mesh_points(a[bl].mesh());
      long const n       = da.extent(1), nw = da.extent(0);
      nda::matrix<double> m(n, n);
      for (long i = 0; i < n; ++i)
        for (long j = i; j < n; ++j) {
          spectral_element(da, i, j, re_a, im_a);
          spectral_element(db, i, j, re_b, im_b);
          re.resize(nw);
          im.resize(nw);
          for (long k = 0; k < nw; ++k) {
            re[k] = re_a[k] - re_b[k];
            im[k] = im_a[k] - im_b[k];
          }
          double const value = nrg::integrate_abs(points, re) + (i == j ? 0.0 : nrg::integrate_abs(points, im));
          m(i, j) = m(j, i) = value;
          result.max        = std::max(result.max, value);
        }
      result.elements.push_back(std::move(m));
    }
    return result;
  }

} // namespace nrgljubljana_interface::dmft
