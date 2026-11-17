// Raw spectral bins of one correlator, as written by NRG Ljubljana (savebins).

#pragma once

#include <h5/h5.hpp>

#include <string>
#include <vector>

namespace nrgljubljana_interface::nrg {

  /// Spectral bins (frequency, weight) of one correlator, in file order.
  /// For a complex spectrum, weight holds the real part and weight_imag the imaginary part of each weight.
  struct spectrum_bins_t {

    /// Frequency of each bin
    std::vector<double> frequency;

    /// Weight of each bin (its real part for a complex spectrum)
    std::vector<double> weight;

    /// Imaginary part of each weight; empty unless is_complex
    std::vector<double> weight_imag;

    /// Whether the spectrum is complex
    bool is_complex = false;

    bool operator==(spectrum_bins_t const &) const = default;
  };

  /// Write spectrum_bins_t to HDF5
  void h5_write(h5::group g, std::string const &name, spectrum_bins_t const &s);

  /// Read spectrum_bins_t from HDF5
  void h5_read(h5::group g, std::string const &name, spectrum_bins_t &s);

  /// Load a binary spectrum file. The layout cannot be detected from the file, so the caller states it:
  /// records of (frequency, weight) for a real spectrum, (frequency, Re weight, Im weight) for a complex one.
  /// A wrong flag is not always caught: a complex file with an even number of records also reads as real.
  /// Throws if the file cannot be read, is not a whole number of records, or contains non-finite values.
  spectrum_bins_t load_spectrum(std::string const &file, bool is_complex);

} // namespace nrgljubljana_interface::nrg
