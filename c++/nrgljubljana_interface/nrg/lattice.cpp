#include "./lattice.hpp"

#include <nrgljubljana/numerics.hpp>

#include <numbers>
#include <stdexcept>

namespace nrgljubljana_interface::nrg {

  std::vector<std::complex<double>> hilbert_moment(dos_table_t const &dos, std::vector<double> const &mesh,
                                                   std::vector<std::complex<double>> const &sigma, double mu, double floor,
                                                   unsigned moment, std::size_t jobs) {
    if (sigma.size() != mesh.size()) throw std::invalid_argument("hilbert_moment: σ has " + std::to_string(sigma.size()) + " values, the mesh " + std::to_string(mesh.size()) + " points");
    NRG::Numerics::Table re{mesh, {}}, im{mesh, {}};
    re.y.reserve(mesh.size());
    im.y.reserve(mesh.size());
    for (auto const &s : sigma) {
      re.y.push_back(s.real());
      im.y.push_back(s.imag());
    }
    auto const l = NRG::Numerics::lattice(NRG::Numerics::Table{dos.energy, dos.density}, re, im, mu, floor, moment, jobs);
    if (l.real.x != mesh or l.spectral.x != mesh) throw std::runtime_error("hilbert_moment: result is not on the mesh");
    std::vector<std::complex<double>> g(mesh.size());
    for (std::size_t k = 0; k < mesh.size(); ++k) g[k] = {l.real.y[k], -std::numbers::pi * l.spectral.y[k]};
    return g;
  }

  double integrate_moment(dos_table_t const &dos, unsigned moment) {
    return NRG::Numerics::integrate_moment(NRG::Numerics::Table{dos.energy, dos.density}, moment);
  }

} // namespace nrgljubljana_interface::nrg
