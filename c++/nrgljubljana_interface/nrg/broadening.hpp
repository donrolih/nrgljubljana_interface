// Broadening of spectral bins, through NRG Ljubljana numerics.

#pragma once

#include "./spectrum.hpp"

#include <vector>

namespace nrgljubljana_interface::nrg {

  /// Parameters of the broadening, as in the broaden tool with -x gamma and omega0_ratio
  struct broadening_t {
    double alpha;               // width of the log-Gaussian
    double temperature;         // T
    double gamma;               // final derivative-of-Fermi pass of width gamma*T
    double omega0_ratio = 1e-9; // omega0 = omega0_ratio*T
  };

  /**
   * Broaden the bins of several twists onto the mesh points and average them over the twists (real part of the
   * bins). The result holds one value per mesh point.
   */
  std::vector<double> broaden(std::vector<spectrum_bins_t const *> const &twists, std::vector<double> const &mesh, broadening_t const &b);

} // namespace nrgljubljana_interface::nrg
