// Reading the raw results of solved twists from their directories.

#pragma once

#include "./options.hpp"
#include "./result.hpp"
#include "../templates/manifest.hpp"

#include <map>
#include <string>

namespace nrgljubljana_interface {

  /// The NRG correlators the options select, with the file of their bins (from the manifest's split of each name)
  std::map<std::string, std::string> selected_spectra(options_t const &r, templates::manifest_t const &m,
                                                       std::map<std::string, char> const &correlator_type);

  /**
   * Read the raw results of one twist from its directory.
   *
   * Always param and td; custom if present (NRG Ljubljana writes it only for templates with singlet operators).
   * With fdm, customfdm and tdfdm. The bins of every given NRG correlator from its file, in the real layout. With energy_flow,
   * annotated.dat. A missing file is an error (std::runtime_error naming the file).
   */
  twist_result_t read_twist(std::string const &dir, double z, std::map<std::string, std::string> const &spectra, bool fdm, bool energy_flow);

} // namespace nrgljubljana_interface
