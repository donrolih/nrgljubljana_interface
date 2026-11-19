// Integrals of functions tabulated on mesh points, through NRG Ljubljana numerics. The function is the Steffen
// interpolant of the values.

#pragma once

#include <vector>

namespace nrgljubljana_interface::nrg {

  /// ∫ f(ω) dω over the mesh
  double integrate(std::vector<double> const &mesh, std::vector<double> const &values);

  /// ∫ |f(ω)| dω over the mesh, split at the roots of the interpolant
  double integrate_abs(std::vector<double> const &mesh, std::vector<double> const &values);

  /// ∫ f(ω) / (exp(ω/T) + 1) dω over the mesh; T = 0 is the step at ω = 0
  double integrate_fermi(std::vector<double> const &mesh, std::vector<double> const &values, double temperature, double epsabs,
                         double epsrel);

} // namespace nrgljubljana_interface::nrg
