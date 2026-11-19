#include "./bethe_lattice.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace nrgljubljana_interface::dmft {

  bethe_lattice_t::bethe_lattice_t(double half_bandwidth, double sigma_floor) : D_(half_bandwidth), floor_(sigma_floor) {
    if (not std::isfinite(D_) or D_ <= 0)
      throw std::invalid_argument("bethe_lattice_t: the half-bandwidth must be positive and finite, got " + detail::number(D_));
    detail::check_sigma_floor(floor_, "bethe_lattice_t");
  }

  detail::scalar_moment_t bethe_lattice_t::scalar() const {
    return [D = D_, floor = floor_](std::string const &, long, std::vector<double> const &mesh, std::vector<std::complex<double>> const &sigma,
                                    double mu, unsigned moment) {
      std::vector<std::complex<double>> g(mesh.size());
      for (std::size_t k = 0; k < mesh.size(); ++k) {
        std::complex<double> const zeta{mesh[k] + mu - sigma[k].real(), std::max(-sigma[k].imag(), floor)};
        std::complex<double> const g0 = 2.0 / (zeta + std::sqrt(zeta - D) * std::sqrt(zeta + D));
        g[k]                          = moment == 0 ? g0 : D * D / 4 * g0 * g0;
      }
      return g;
    };
  }

  std::pair<g_w_t, g_w_t> bethe_lattice_t::moments(g_w_t const &Sigma, double mu) const {
    auto m = detail::lattice_moments(Sigma, mu, 1, scalar(), true, "bethe_lattice_t");
    return {std::move(m[0]), std::move(m[1])};
  }

  g_w_t bethe_lattice_t::local_gf(g_w_t const &Sigma, double mu) const {
    return std::move(detail::lattice_moments(Sigma, mu, 0, scalar(), true, "bethe_lattice_t")[0]);
  }

  block_matrices_t bethe_lattice_t::onsite_energy(triqs::gfs::gf_struct_t const &structure) const {
    block_matrices_t result;
    for (auto const &[name, n] : structure) {
      nda::matrix<std::complex<double>> m(n, n);
      m = 0;
      result.push_back(std::move(m));
    }
    return result;
  }

} // namespace nrgljubljana_interface::dmft
