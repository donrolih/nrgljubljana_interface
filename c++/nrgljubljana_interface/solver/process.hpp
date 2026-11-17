// From the raw results of the twists to Green's functions on the solver mesh.

#pragma once

#include "./result.hpp"
#include "../templates/manifest.hpp"

#include <triqs/mesh.hpp>

namespace nrgljubljana_interface {

  /**
   * Fill the averaged and processed parts of a result from its twists.
   *
   * Averages: expv and tdfdm, the mean over the twists of customfdm and tdfdm (without T).
   *
   * Processed (if the request broadens): for each requested logical correlator a Green's function on mesh with the
   * structure from the manifest. Each listed element is the broadened spectral density A of its NRG correlator,
   * averaged over the twists, as Im = -pi A for fermionic types (d, q) and Im = A for bosonic types (s, t, o), with
   * the real part by Kramers-Kronig; unlisted elements are zero. With G, F_l, F_r and I also "SigmaHartree", the
   * twist average of the [hartree] expectation values (0 where none), and "Sigma" = SigmaHartree + I - F_l G^-1 F_r.
   *
   * Broadening: alpha, T and gamma of the solve parameters, omega0 = 1e-9 T, as the broaden tool with -x gamma.
   */
  void process(solve_result_t &r, templates::manifest_t const &m, triqs::mesh::refreq_log const &mesh);

} // namespace nrgljubljana_interface
