#include "./collect.hpp"

#include "../nrg/output.hpp"
#include "../nrg/spectrum.hpp"

#include <filesystem>
#include <stdexcept>

namespace nrgljubljana_interface {

  std::map<std::string, char> selected_spectra(options_t const &r, templates::manifest_t const &m,
                                                std::map<std::string, char> const &correlator_type) {
    std::map<std::string, char> result;
    for (auto const &name : selected_correlators(r, m))
      for (auto const &entry : m.correlators.at(name)) result[entry.nrg_name] = correlator_type.at(entry.nrg_name);
    return result;
  }

  twist_result_t read_twist(std::string const &dir, double z, std::map<std::string, char> const &spectra, bool fdm, bool energy_flow) {
    auto const path   = [&](std::string const &file) { return (std::filesystem::path(dir) / file).string(); };
    auto const exists = [&](std::string const &file) { return std::filesystem::exists(path(file)); };

    twist_result_t t;
    t.z     = z;
    t.param = nrg::read_text_file(path("param"));
    t.td    = nrg::parse_step_table(nrg::read_text_file(path("td")), path("td"));
    if (exists("custom")) t.custom = nrg::parse_step_table(nrg::read_text_file(path("custom")), path("custom"));
    if (fdm) {
      t.customfdm = nrg::parse_final_values(nrg::read_text_file(path("customfdm")), path("customfdm"));
      t.tdfdm     = nrg::parse_final_values(nrg::read_text_file(path("tdfdm")), path("tdfdm"));
    }
    for (auto const &[corr, type] : spectra) {
      auto const file = path(nrg::spectrum_file(corr, type));
      if (not std::filesystem::exists(file)) throw std::runtime_error(file + ": missing; NRG Ljubljana did not compute " + corr);
      t.spectra[corr] = nrg::load_spectrum(file, false);
    }
    if (energy_flow) t.energy_flow = nrg::parse_energy_flow(nrg::read_text_file(path("annotated.dat")), path("annotated.dat"));
    return t;
  }

} // namespace nrgljubljana_interface
