#include "./options.hpp"

#include <algorithm>
#include <stdexcept>

namespace nrgljubljana_interface {

  void validate_options(options_t const &r, templates::manifest_t const &m) {
    auto const fail = [](std::string const &what) { throw std::invalid_argument("Invalid solve options: " + what); };

    if (r.correlators) {
      if (not r.fdm and not r.correlators->empty()) fail("correlators require fdm=true");
      auto const &names = *r.correlators;
      for (std::size_t k = 0; k < names.size(); ++k) {
        if (not m.correlators.contains(names[k])) fail("unknown correlator \"" + names[k] + "\" (not in the template's manifest)");
        if (std::find(names.begin(), names.begin() + static_cast<std::ptrdiff_t>(k), names[k]) != names.begin() + static_cast<std::ptrdiff_t>(k))
          fail("correlator \"" + names[k] + "\" given twice");
      }
    }
    if (not r.broaden and not r.keep_bins) fail("broaden=false requires keep_bins=true");
    if (r.energies < 0) fail("energies must not be negative");
    for (std::size_t k = 0; k < r.z.size(); ++k) {
      auto const z = r.z[k];
      if (not(z > 0.0 and z <= 1.0)) fail("z = " + std::to_string(z) + " is outside (0, 1]");
      if (std::find(r.z.begin(), r.z.begin() + static_cast<std::ptrdiff_t>(k), z) != r.z.begin() + static_cast<std::ptrdiff_t>(k))
        fail("z = " + std::to_string(z) + " given twice");
    }
  }

  std::vector<std::string> selected_correlators(options_t const &r, templates::manifest_t const &m) {
    if (r.correlators) return *r.correlators;
    std::vector<std::string> all;
    if (r.fdm)
      for (auto const &[name, entries] : m.correlators) all.push_back(name);
    return all;
  }

  void h5_write(h5::group g, std::string const &name, options_t const &r) {
    auto grp = g.create_group(name);
    h5_write(grp, "fdm", r.fdm);
    h5_write(grp, "correlators", r.correlators);
    h5_write(grp, "broaden", r.broaden);
    h5_write(grp, "keep_bins", r.keep_bins);
    h5_write(grp, "energies", r.energies);
    h5_write(grp, "z", r.z);
    h5_write(grp, "keep_workdir", r.keep_workdir);
  }

  void h5_read(h5::group g, std::string const &name, options_t &r) {
    auto grp = g.open_group(name);
    h5_read(grp, "fdm", r.fdm);
    r.correlators.reset();
    h5_read(grp, "correlators", r.correlators);
    h5_read(grp, "broaden", r.broaden);
    h5_read(grp, "keep_bins", r.keep_bins);
    h5_read(grp, "energies", r.energies);
    h5_read(grp, "z", r.z);
    h5_read(grp, "keep_workdir", r.keep_workdir);
  }

} // namespace nrgljubljana_interface
