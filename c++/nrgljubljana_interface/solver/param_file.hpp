// The NRG Ljubljana parameter file of a solve, and how a request sets it.

#pragma once

#include "./request.hpp"
#include "../params.hpp"
#include "../templates/manifest.hpp"

#include <map>
#include <string>
#include <vector>

namespace nrgljubljana_interface {

  /// The parameter file ([extra] model parameters and [param]) for the given parameters
  std::string render_param_file(constr_params_t const &cp, solve_params_t const &sp, nrg_params_t const &np);

  /**
   * Set the parameters a request decides.
   *
   * Spectrum lists (specs, specd, spect, specq, specot): the NRG correlators of the requested logical correlators,
   * by type; savebins on if any are requested. fdm=false: fdm, fdmexpv, dmnrg and cfs off. fdm=true: fdm and fdmexpv
   * on. energies > 0: dumpannotated=energies; otherwise dumpannotated is left as is.
   */
  void apply_request(solve_request_t const &r, templates::manifest_t const &m, std::map<std::string, char> const &correlator_type,
                     constr_params_t &cp, nrg_params_t &np);

  /// The twist values a request solves: its z, or i/Nz for i = 1..Nz
  std::vector<double> twist_values(solve_request_t const &r, int Nz);

  /// Whether z is the full set i/Nz, i = 1..Nz (in any order)
  bool complete_z_average(std::vector<double> const &z, int Nz);

  /**
   * Check the deprecated solve_params_t::method against a request.
   *
   * "" (not set) and "fdm" are accepted; "fdm" with a request without fdm is an error, as is any other value.
   * Returns true if method was set, so that the caller can warn that it is deprecated.
   */
  bool check_method(solve_params_t const &sp, solve_request_t const &r);

} // namespace nrgljubljana_interface
