// What a solve computes and returns.

#pragma once

#include "../templates/manifest.hpp"

#include <h5/h5.hpp>

#include <optional>
#include <string>
#include <vector>

namespace nrgljubljana_interface {

  /// A solve request: what NRG Ljubljana computes and what the solver returns
  struct solve_request_t {

    /// Run the FDM sweep after the forward sweep (FDM expectation values, thermodynamics and spectra)
    bool fdm = true;

    /// Logical correlators to compute (G, F_l, F_r, I, chi_SS, ...: names from the manifest).
    /// Unset: all the template provides if fdm, none otherwise. Empty: none.
    std::optional<std::vector<std::string>> correlators = {};

    /// Broaden the requested spectra in the solve; without it the result holds only the bins
    bool broaden = true;

    /// Keep the raw spectral bins in the result; required if not broaden
    bool keep_bins = true;

    /// Lowest levels per NRG step returned in the energy flow; 0 for none
    long energies = 0;

    /// Twist values in (0, 1] to solve; empty for i/Nz, i = 1..Nz
    std::vector<double> z = {};

    /// Keep the working directory after the solve
    bool keep_workdir = false;

    bool operator==(solve_request_t const &) const = default;
  };

  /**
   * Check a request against a template's manifest.
   *
   * Throws std::invalid_argument for correlators requested without fdm, broaden=false without keep_bins, an unknown
   * correlator or one given twice, energies < 0, and a z value outside (0, 1] or given twice.
   */
  void validate_request(solve_request_t const &r, templates::manifest_t const &m);

  /// The logical correlators a request computes: as listed, or by default all of the manifest's if fdm, else none
  std::vector<std::string> requested_correlators(solve_request_t const &r, templates::manifest_t const &m);

  /// Write solve_request_t to HDF5
  void h5_write(h5::group g, std::string const &name, solve_request_t const &r);

  /// Read solve_request_t from HDF5
  void h5_read(h5::group g, std::string const &name, solve_request_t &r);

} // namespace nrgljubljana_interface
