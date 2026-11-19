// Kramers-Kronig transforms, through NRG Ljubljana numerics.

#pragma once

#include <vector>

namespace nrgljubljana_interface::nrg {

  /// Kramers-Kronig transform of values on a mesh symmetric about zero: Im part to Re part, as the kk tool
  std::vector<double> kramers_kronig(std::vector<double> const &mesh, std::vector<double> const &values);

} // namespace nrgljubljana_interface::nrg
