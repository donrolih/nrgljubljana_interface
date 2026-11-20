#include "./callable_lattice.hpp"

#include <cmath>
#include <stdexcept>

namespace nrgljubljana_interface::dmft {

  namespace {

    void check_result(g_w_t const &G, g_w_t const &Sigma, std::string const &what) {
      detail::check_same_structure(G, Sigma, "callable_lattice_t", what, "Sigma");
      detail::check_finite(G, "callable_lattice_t", what);
    }

  } // namespace

  callable_lattice_t::callable_lattice_t(moments_t moments, block_matrices_t onsite_energy, std::optional<local_gf_t> local_gf, std::string description)
     : moments_(std::move(moments)), onsite_(std::move(onsite_energy)), local_gf_(std::move(local_gf)), description_(std::move(description)) {
    if (not moments_) throw std::invalid_argument("callable_lattice_t: no moments function");
    if (local_gf_ and not *local_gf_) throw std::invalid_argument("callable_lattice_t: an empty local_gf function");
    for (std::size_t b = 0; b < onsite_.size(); ++b) {
      auto const &m = onsite_[b];
      if (m.extent(0) != m.extent(1))
        throw std::invalid_argument("callable_lattice_t: the on-site energy of block " + std::to_string(b) + " is not square");
      for (long i = 0; i < m.extent(0); ++i)
        for (long j = 0; j < m.extent(1); ++j)
          if (not std::isfinite(m(i, j).real()) or not std::isfinite(m(i, j).imag()))
            throw std::invalid_argument("callable_lattice_t: the on-site energy of block " + std::to_string(b) + " is not finite");
    }
  }

  std::pair<g_w_t, g_w_t> callable_lattice_t::moments(g_w_t const &Sigma, double mu) const {
    auto result = moments_(Sigma, mu);
    check_result(result.first, Sigma, "G0");
    check_result(result.second, Sigma, "G1");
    return result;
  }

  g_w_t callable_lattice_t::local_gf(g_w_t const &Sigma, double mu) const {
    if (not local_gf_) return moments(Sigma, mu).first;
    auto G = (*local_gf_)(Sigma, mu);
    check_result(G, Sigma, "G0");
    return G;
  }

  block_matrices_t callable_lattice_t::onsite_energy(triqs::gfs::gf_struct_t const &structure) const {
    if (structure.size() != onsite_.size())
      throw std::invalid_argument("callable_lattice_t: the on-site energy has " + std::to_string(onsite_.size()) + " blocks, the structure "
                                  + std::to_string(structure.size()));
    for (std::size_t b = 0; b < structure.size(); ++b) {
      auto const &[name, n] = structure[b];
      if (onsite_[b].extent(0) != n)
        throw std::invalid_argument("callable_lattice_t: the on-site energy of block " + name + " is " + std::to_string(onsite_[b].extent(0)) + "x"
                                    + std::to_string(onsite_[b].extent(0)) + ", the block has size " + std::to_string(n));
    }
    return onsite_;
  }

} // namespace nrgljubljana_interface::dmft
