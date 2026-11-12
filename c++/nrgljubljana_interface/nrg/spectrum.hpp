// Raw spectral bins of one correlator, as written by NRG Ljubljana (savebins).

#pragma once

#include <string>
#include <vector>

namespace nrgljubljana_interface::nrg {

  /// Spectral bins (frequency, weight) of one correlator, in file order.
  /// For a complex spectrum, weight holds the real part and weight_imag the imaginary part of each weight.
  struct spectrum_bins_t {
    std::vector<double> frequency;
    std::vector<double> weight;
    std::vector<double> weight_imag; // empty unless is_complex
    bool is_complex = false;
  };

  /// Load a binary spectrum file. The layout cannot be detected from the file, so the caller states it:
  /// records of (frequency, weight) for a real spectrum, (frequency, Re weight, Im weight) for a complex one.
  /// A wrong flag is not always caught: a complex file with an even number of records also reads as real.
  /// Throws if the file cannot be read, is not a whole number of records, or contains non-finite values.
  spectrum_bins_t load_spectrum(std::string const &file, bool is_complex);

} // namespace nrgljubljana_interface::nrg
