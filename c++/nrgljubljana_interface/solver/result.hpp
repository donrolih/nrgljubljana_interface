// The result of a solve: raw data of each twist, their averages, and the processed Green's functions.

#pragma once

#include "./options.hpp"
#include "../nrg/chain.hpp"
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

    /// Twist value
    double z = 1;

    /// The rendered NRG Ljubljana parameter file
    std::string param;

    /// Spectral bins of each NRG correlator
    std::map<std::string, nrg::spectrum_bins_t> spectra;

    /// Expectation values per NRG step (custom)
    step_table_t custom;

    /// Thermodynamics per NRG step (td)
    step_table_t td;

    /// FDM expectation values (customfdm)
    std::map<std::string, double> customfdm;

    /// FDM thermodynamics (tdfdm)
    std::map<std::string, double> tdfdm;

    /// Lowest levels of each NRG step (empty unless options.energies > 0)
    energy_flow_t energy_flow;
    /// The Wilson chain coefficients of the twist (in-process route; empty on the script route)
    nrg::chain_coefficients_t chain = {};
  };

  /// The result of a solve: what was solved, the raw data of each twist, their averages and the Green's functions
  struct result_t {

    /// Options of the solve
    options_t options;

    /// Solve parameters, as given
    solve_params_t solve_params;

    /// NRG parameters, as given (the values used for each twist are in its parameter file)
    nrg_params_t nrg_params;

    /// Model (template)
    std::string model;

    /// Symmetry type (template)
    std::string symtype;

    /// The manifest text of the template
    std::string manifest;

    /// Spectrum list of each NRG correlator: s, t, o, d or q
    std::map<std::string, std::string> correlator_type;

    /// NRG Ljubljana source revision
    std::string nrg_revision;

    /// The solver mesh: of processed, and the default of broaden
    triqs::mesh::refreq_log mesh;

    /// Twist values solved
    std::vector<double> z;

    /// Whether z is the full set i/Nz
    bool complete_z_average = false;

    /// Raw data of each twist
    std::vector<twist_result_t> twists;

    /// The points (block, frequency) at which an eigenvalue of Γ was raised to the floor, and the largest amount
    /// (in-process route; 0 on the script route)
    long Gamma_projected_points = 0;
    double Gamma_largest_violation = 0;

    /// FDM expectation values, averaged over the twists
    std::map<std::string, double> expv;

    /// FDM thermodynamics, averaged over the twists
    std::map<std::string, double> tdfdm;

    /// Green's function of each broadened correlator (and Sigma, SigmaHartree) on the solver mesh
    std::map<std::string, g_w_t> processed;

    /// A copy without the raw spectral bins: smaller in HDF5, but cannot be broadened again
    [[nodiscard]] result_t without_bins() const;

    /// HDF5 format tag
    static std::string hdf5_format() { return "NRGLJUBLJANA_INTERFACE_Result"; }
  };

  /// Write result_t to HDF5, with the raw spectral bins
  void h5_write(h5::group g, std::string const &name, result_t const &r);

  /// Write result_t to HDF5, with or without the raw spectral bins
  void h5_write(h5::group g, std::string const &name, result_t const &r, bool with_bins);

  /// Read result_t from HDF5; a result written without bins reads with empty spectra
  void h5_read(h5::group g, std::string const &name, result_t &r);

} // namespace nrgljubljana_interface
