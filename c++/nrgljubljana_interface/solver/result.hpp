// The result of a solve: raw data of each twist, their averages, and the processed Green's functions.

#pragma once

#include "./request.hpp"
#include "../nrg/output.hpp"
#include "../nrg/spectrum.hpp"
#include "../params.hpp"
#include "../types.hpp"

#include <h5/h5.hpp>
#include <nda/nda.hpp>

#include <map>
#include <string>
#include <vector>

namespace nrgljubljana_interface {

  using nrg::energy_flow_t;
  using nrg::energy_level_t;
  using nrg::step_table_t;

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
