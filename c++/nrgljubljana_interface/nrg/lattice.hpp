// Hilbert transforms of a tabulated density of states, through NRG Ljubljana numerics.

#pragma once

#include <complex>
#include <cstddef>
#include <vector>

namespace nrgljubljana_interface::nrg {

  /// A tabulated density of states: strictly increasing energies, non-negative densities
  struct dos_table_t {
    std::vector<double> energy, density;
  };

  /**
   * G_n(ω) = ∫ εⁿ ρ(ε) / (z(ω) − ε) dε at the mesh points, with z = ω + μ − Re σ(ω) + i max(−Im σ(ω), floor) and the
   * Steffen interpolant of ρ. σ holds one value per mesh point; the result too.
   */
  std::vector<std::complex<double>> hilbert_moment(dos_table_t const &dos, std::vector<double> const &mesh,
                                                   std::vector<std::complex<double>> const &sigma, double mu, double floor,
                                                   unsigned moment, std::size_t jobs);

  /// ∫ εⁿ ρ(ε) dε of the Steffen interpolant
  double integrate_moment(dos_table_t const &dos, unsigned moment);

} // namespace nrgljubljana_interface::nrg
