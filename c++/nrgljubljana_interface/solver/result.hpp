// The result of a solve: raw data of each twist, their averages, and the processed Green's functions.

#pragma once

#include "./request.hpp"
#include "../nrg/spectrum.hpp"
#include "../params.hpp"
#include "../types.hpp"

#include <h5/h5.hpp>
#include <nda/nda.hpp>

#include <map>
#include <string>
#include <vector>

namespace nrgljubljana_interface {

  /// A table with one row per NRG step (custom, td): column names and values
  struct step_table_t {
    std::vector<std::string> columns;
    nda::array<double, 2> values; // (step, column)
  };

  /// A group of degenerate levels of one NRG step, as in annotated.dat
  struct energy_level_t {
    double energy = 0;
    std::vector<std::vector<int>> subspaces; // quantum numbers of each invariant subspace in the group
    long degeneracy = 0;                     // total number of states, including multiplicity

    bool operator==(energy_level_t const &) const = default;
  };

  /// Lowest levels of each NRG step
  using energy_flow_t = std::vector<std::vector<energy_level_t>>;

  /// Raw data of one twist
  struct twist_result_t {
    double z = 1;
    std::string param;                                     // the rendered parameter file
    std::map<std::string, nrg::spectrum_bins_t> spectra;   // NRG correlator -> bins
    step_table_t custom, td;                               // per step
    std::map<std::string, double> customfdm, tdfdm;        // FDM values
    energy_flow_t energy_flow;
  };

  struct solve_result_t {
    // What was solved
    solve_request_t request;
    solve_params_t solve_params;
    nrg_params_t nrg_params;
    std::string model, symtype;
    std::string manifest;                       // the template's manifest text
    std::map<std::string, char> correlator_type; // NRG correlator -> spectrum list (s, t, o, d, q)
    std::string nrg_revision;                   // NRG Ljubljana source revision

    // The twists
    std::vector<double> z;
    bool complete_z_average = false; // z is the full set i/Nz
    std::vector<twist_result_t> twists;

    // Averaged over the twists
    std::map<std::string, double> expv, tdfdm;

    // Logical correlator (and Sigma) -> Green's function on the solver mesh
    std::map<std::string, g_w_t> processed;
  };

  /// Write solve_result_t to HDF5, with the raw spectral bins
  void h5_write(h5::group g, std::string const &name, solve_result_t const &r);

  /// Write solve_result_t to HDF5, with or without the raw spectral bins
  void h5_write(h5::group g, std::string const &name, solve_result_t const &r, bool with_bins);

  /// Read solve_result_t from HDF5; a result written without bins reads with empty spectra
  void h5_read(h5::group g, std::string const &name, solve_result_t &r);

} // namespace nrgljubljana_interface
