// What a solve computes and returns.

#pragma once

#include "../templates/manifest.hpp"

#include <h5/h5.hpp>

#include <optional>
#include <string>
#include <vector>

namespace nrgljubljana_interface {

  /// Options of a solve: what NRG Ljubljana computes and what the solver returns
  struct options_t {

    /// Run the FDM sweep after the forward sweep (FDM expectation values, thermodynamics and spectra)
    bool fdm = true;

    /// Correlators to compute (G, F_l, F_r, I, chi_SS, ...: names from the manifest).
    /// Unset: all the template provides if fdm, none otherwise. Empty: none.
    std::optional<std::vector<std::string>> correlators = {};

    /// Broaden the selected spectra in the solve; without it the result holds only the bins
    bool broaden = true;

    /// Keep the raw spectral bins in the result; required if not broaden
    bool keep_bins = true;

    /// Lowest levels per NRG step returned in the energy flow; 0 for none
    long energies = 0;

    /// Twist values in (0, 1] to solve; empty for i/Nz, i = 1..Nz
    std::vector<double> z = {};

    /// Keep the working directory after the solve
    bool keep_workdir = false;

    bool operator==(options_t const &) const = default;
  };

  /**
   * Check options against a template's manifest.
   *
   * Throws std::invalid_argument for correlators selected without fdm, broaden=false without keep_bins, an unknown
   * correlator or one given twice, energies < 0, and a z value outside (0, 1] or given twice.
   */
  void validate_options(options_t const &r, templates::manifest_t const &m);

  /// The correlators the options select: as listed, or by default all of the manifest's if fdm, else none
  std::vector<std::string> selected_correlators(options_t const &r, templates::manifest_t const &m);

  /// Write options_t to HDF5
  void h5_write(h5::group g, std::string const &name, options_t const &r);

  /// Read options_t from HDF5
  void h5_read(h5::group g, std::string const &name, options_t &r);

} // namespace nrgljubljana_interface
