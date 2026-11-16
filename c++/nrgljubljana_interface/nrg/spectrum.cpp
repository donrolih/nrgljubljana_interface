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

  void h5_write(h5::group g, std::string const &name, spectrum_bins_t const &s) {
    auto grp = g.create_group(name);
    h5_write(grp, "frequency", s.frequency);
    h5_write(grp, "weight", s.weight);
    h5_write(grp, "is_complex", s.is_complex);
    if (s.is_complex) h5_write(grp, "weight_imag", s.weight_imag);
  }

  void h5_read(h5::group g, std::string const &name, spectrum_bins_t &s) {
    auto grp = g.open_group(name);
    h5_read(grp, "frequency", s.frequency);
    h5_read(grp, "weight", s.weight);
    h5_read(grp, "is_complex", s.is_complex);
    s.weight_imag.clear();
    if (s.is_complex) h5_read(grp, "weight_imag", s.weight_imag);
  }

} // namespace nrgljubljana_interface::nrg
