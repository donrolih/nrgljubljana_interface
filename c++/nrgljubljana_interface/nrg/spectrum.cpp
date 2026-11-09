#include "./spectrum.hpp"

#include <nrgljubljana/numerics.hpp>

#include <cstddef>
#include <stdexcept>

namespace nrgljubljana_interface::nrg {

  spectrum_bins_t load_spectrum(std::string const &file, bool is_complex) {
    using NRG::Numerics::SpectrumLayout;
    spectrum_bins_t s;
    s.is_complex = is_complex;

    // NRG Ljubljana reads one weight column per call: the real part here, the imaginary part below
    auto const re = NRG::Numerics::load_spectrum_bins(file, is_complex ? SpectrumLayout::complex_real : SpectrumLayout::real_pairs);
    s.frequency.reserve(re.size());
    s.weight.reserve(re.size());
    for (auto const &b : re) {
      s.frequency.push_back(b.frequency);
      s.weight.push_back(b.weight);
    }

    if (is_complex) {
      auto const im = NRG::Numerics::load_spectrum_bins(file, SpectrumLayout::complex_imaginary);
      // Same file, so the same frequencies, unless it changed between the two reads
      if (im.size() != re.size()) throw std::runtime_error("Binary spectrum changed while it was read: " + file);
      s.weight_imag.reserve(im.size());
      for (std::size_t k = 0; k < im.size(); ++k) {
        if (im[k].frequency != re[k].frequency) throw std::runtime_error("Binary spectrum changed while it was read: " + file);
        s.weight_imag.push_back(im[k].weight);
      }
    }
    return s;
  }

} // namespace nrgljubljana_interface::nrg
