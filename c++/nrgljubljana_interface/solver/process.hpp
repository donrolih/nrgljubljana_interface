// From the raw results of the twists to Green's functions on the solver mesh.

#pragma once

#include "./result.hpp"
#include "../templates/manifest.hpp"

#include <triqs/mesh.hpp>

#include <optional>
#include <string>

namespace nrgljubljana_interface {

  /// Parameters of a broadening of stored bins; T is that of the solve
  struct broadening_params_t {

    /// Width of the log-Gaussian kernel; unset: alpha of the solve
    std::optional<double> alpha = {};

    /// Width of the final derivative-of-Fermi pass, in units of T; unset: gamma of the solve
    std::optional<double> gamma = {};

    /// omega0 = omega0_ratio * T
    double omega0_ratio = 1e-9;

    bool operator==(broadening_params_t const &) const = default;
  };

  /**
   * Fill the averaged and processed parts of a result from its twists, on the mesh of the result.
   *
   * Averages: expv and tdfdm, the mean over the twists of customfdm and tdfdm (without T).
   *
   * Processed (if the options broaden): each selected logical correlator, as broaden(r, name) with the parameters
   * of the solve. With G, F_l, F_r and I also "SigmaHartree", the twist average of the [hartree] expectation values
   * (0 where none, constant in frequency), and "Sigma" = SigmaHartree + I - F_l G^-1 F_r.
   */
  void process(result_t &r, templates::manifest_t const &m);

  /**
   * Broaden the stored bins of a logical correlator of a result into a Green's function, without running NRG
   * Ljubljana.
   *
   * The Green's function has the structure from the manifest. Each listed element is the broadened spectral density
   * A of its NRG correlator, averaged over the twists of the result, as Im = -pi A for fermionic types (d, q) and
   * Im = A for bosonic types (s, t, o), with the real part by Kramers-Kronig; unlisted elements are zero.
   *
   * The mesh defaults to the mesh of the result. Throws std::invalid_argument for a name not in the manifest (Sigma
   * and SigmaHartree are not correlators) and std::runtime_error if the result holds no bins of the correlator.
   */
  g_w_t broaden(result_t const &r, std::string const &correlator, broadening_params_t const &b = {},
                std::optional<triqs::mesh::refreq_log> const &mesh = {});

} // namespace nrgljubljana_interface
