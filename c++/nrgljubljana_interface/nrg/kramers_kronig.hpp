// Kramers-Kronig transforms, through NRG Ljubljana numerics.

#pragma once

#include <cstddef>
#include <vector>

namespace nrgljubljana_interface::nrg {

  /// Kramers-Kronig transform of values on a mesh symmetric about zero: Im part to Re part, as the kk tool, with
  /// the given number of workers (the result does not depend on it)
  std::vector<double> kramers_kronig(std::vector<double> const &mesh, std::vector<double> const &values, std::size_t jobs = 1);

} // namespace nrgljubljana_interface::nrg
