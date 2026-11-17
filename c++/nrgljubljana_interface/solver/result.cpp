#include "./result.hpp"

#include <nda/h5.hpp>

#include <stdexcept>

namespace nrgljubljana_interface {

  namespace {

    void write_table(h5::group g, std::string const &name, step_table_t const &t) {
      auto grp = g.create_group(name);
      h5_write(grp, "columns", t.columns);
      h5_write(grp, "values", t.values);
    }

    void read_table(h5::group g, std::string const &name, step_table_t &t) {
      auto grp = g.open_group(name);
      h5_read(grp, "columns", t.columns);
      h5_read(grp, "values", t.values);
    }

    // Stored flat: levels per step, then per level its energy, degeneracy and subspaces, then per subspace its
    // quantum numbers
    void write_energy_flow(h5::group g, std::string const &name, energy_flow_t const &flow) {
      std::vector<long> step_levels, level_subspaces, subspace_size, degeneracy;
      std::vector<double> energy;
      std::vector<int> quantum_numbers;
      for (auto const &step : flow) {
        step_levels.push_back(static_cast<long>(step.size()));
        for (auto const &level : step) {
          energy.push_back(level.energy);
          degeneracy.push_back(level.degeneracy);
          level_subspaces.push_back(static_cast<long>(level.subspaces.size()));
          for (auto const &qn : level.subspaces) {
            subspace_size.push_back(static_cast<long>(qn.size()));
            quantum_numbers.insert(quantum_numbers.end(), qn.begin(), qn.end());
          }
        }
      }
      auto grp = g.create_group(name);
      h5_write(grp, "step_levels", step_levels);
      h5_write(grp, "energy", energy);
      h5_write(grp, "degeneracy", degeneracy);
      h5_write(grp, "level_subspaces", level_subspaces);
      h5_write(grp, "subspace_size", subspace_size);
      h5_write(grp, "quantum_numbers", quantum_numbers);
    }

    void read_energy_flow(h5::group g, std::string const &name, energy_flow_t &flow) {
      std::vector<long> step_levels, level_subspaces, subspace_size, degeneracy;
      std::vector<double> energy;
      std::vector<int> quantum_numbers;
      auto grp = g.open_group(name);
      h5_read(grp, "step_levels", step_levels);
      h5_read(grp, "energy", energy);
      h5_read(grp, "degeneracy", degeneracy);
      h5_read(grp, "level_subspaces", level_subspaces);
      h5_read(grp, "subspace_size", subspace_size);
      h5_read(grp, "quantum_numbers", quantum_numbers);

      auto const corrupt = [&] { throw std::runtime_error("Corrupt energy flow in HDF5 group " + name); };
      if (energy.size() != degeneracy.size() or energy.size() != level_subspaces.size()) corrupt();
      flow.clear();
      std::size_t level = 0, subspace = 0, qn = 0;
      for (auto const n_levels : step_levels) {
        auto &step = flow.emplace_back();
        for (long l = 0; l < n_levels; ++l, ++level) {
          if (level >= energy.size()) corrupt();
          energy_level_t e{energy[level], {}, degeneracy[level]};
          for (long s = 0; s < level_subspaces[level]; ++s, ++subspace) {
            if (subspace >= subspace_size.size() or qn + subspace_size[subspace] > quantum_numbers.size()) corrupt();
            auto const first = quantum_numbers.begin() + static_cast<std::ptrdiff_t>(qn);
            e.subspaces.emplace_back(first, first + subspace_size[subspace]);
            qn += subspace_size[subspace];
          }
          step.push_back(std::move(e));
        }
      }
      if (level != energy.size() or subspace != subspace_size.size() or qn != quantum_numbers.size()) corrupt();
    }

    void write_twist(h5::group g, std::string const &name, twist_result_t const &t, bool with_bins) {
      auto grp = g.create_group(name);
      h5_write(grp, "z", t.z);
      h5_write(grp, "param", t.param);
      if (with_bins) h5_write(grp, "spectra", t.spectra);
      write_table(grp, "custom", t.custom);
      write_table(grp, "td", t.td);
      h5_write(grp, "customfdm", t.customfdm);
      h5_write(grp, "tdfdm", t.tdfdm);
      write_energy_flow(grp, "energy_flow", t.energy_flow);
    }

    void read_twist(h5::group g, std::string const &name, twist_result_t &t) {
      auto grp = g.open_group(name);
      h5_read(grp, "z", t.z);
      h5_read(grp, "param", t.param);
      t.spectra.clear();
      if (grp.has_key("spectra")) h5_read(grp, "spectra", t.spectra);
      read_table(grp, "custom", t.custom);
      read_table(grp, "td", t.td);
      h5_read(grp, "customfdm", t.customfdm);
      h5_read(grp, "tdfdm", t.tdfdm);
      read_energy_flow(grp, "energy_flow", t.energy_flow);
    }

  } // namespace

  result_t result_t::without_bins() const {
    auto r = *this;
    for (auto &t : r.twists) t.spectra.clear();
    return r;
  }

  void h5_write(h5::group g, std::string const &name, result_t const &r) { h5_write(g, name, r, true); }

  void h5_write(h5::group g, std::string const &name, result_t const &r, bool with_bins) {
    auto grp = g.create_group(name);
    h5::write_hdf5_format(grp, r);
    h5_write(grp, "options", r.options);
    h5_write(grp, "solve_params", r.solve_params);
    h5_write(grp, "nrg_params", r.nrg_params);
    h5_write(grp, "model", r.model);
    h5_write(grp, "symtype", r.symtype);
    h5_write(grp, "manifest", r.manifest);
    h5_write(grp, "correlator_type", r.correlator_type);
    h5_write(grp, "nrg_revision", r.nrg_revision);
    h5_write(grp, "mesh", r.mesh);
    h5_write(grp, "z", r.z);
    h5_write(grp, "complete_z_average", r.complete_z_average);
    auto twists = grp.create_group("twists");
    for (std::size_t k = 0; k < r.twists.size(); ++k) write_twist(twists, std::to_string(k), r.twists[k], with_bins);
    h5_write(grp, "expv", r.expv);
    h5_write(grp, "tdfdm", r.tdfdm);
    h5_write(grp, "processed", r.processed);
  }

  void h5_read(h5::group g, std::string const &name, result_t &r) {
    auto grp = g.open_group(name);
    h5_read(grp, "options", r.options);
    h5_read(grp, "solve_params", r.solve_params);
    h5_read(grp, "nrg_params", r.nrg_params);
    h5_read(grp, "model", r.model);
    h5_read(grp, "symtype", r.symtype);
    h5_read(grp, "manifest", r.manifest);
    r.correlator_type.clear();
    h5_read(grp, "correlator_type", r.correlator_type);
    h5_read(grp, "nrg_revision", r.nrg_revision);
    h5_read(grp, "mesh", r.mesh);
    h5_read(grp, "z", r.z);
    h5_read(grp, "complete_z_average", r.complete_z_average);
    auto twists = grp.open_group("twists");
    r.twists.clear();
    for (std::size_t k = 0; twists.has_key(std::to_string(k)); ++k) read_twist(twists, std::to_string(k), r.twists.emplace_back());
    h5_read(grp, "expv", r.expv);
    h5_read(grp, "tdfdm", r.tdfdm);
    r.processed.clear();
    h5_read(grp, "processed", r.processed);
  }

} // namespace nrgljubljana_interface
