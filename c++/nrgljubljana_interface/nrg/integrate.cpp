#include "./integrate.hpp"

#include <nrgljubljana/numerics.hpp>

namespace nrgljubljana_interface::nrg {

  double integrate(std::vector<double> const &mesh, std::vector<double> const &values) {
    return NRG::Numerics::integrate(NRG::Numerics::Table{mesh, values});
  }

  double integrate_abs(std::vector<double> const &mesh, std::vector<double> const &values) {
    return NRG::Numerics::integrate_abs(NRG::Numerics::Table{mesh, values});
  }

  double integrate_fermi(std::vector<double> const &mesh, std::vector<double> const &values, double temperature, double epsabs,
                         double epsrel) {
    return NRG::Numerics::integrate_fermi(NRG::Numerics::Table{mesh, values}, temperature, epsabs, epsrel);
  }

} // namespace nrgljubljana_interface::nrg
