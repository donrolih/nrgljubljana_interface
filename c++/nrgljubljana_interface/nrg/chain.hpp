// The in-process route of a solve through NRG Ljubljana's workflow library: the Wilson chains of Γ (mixchain), the
// initialization from their coefficients, and the solve of one twist.

#pragma once

#include <h5/h5.hpp>

#include <map>
#include <string>
#include <vector>

namespace nrgljubljana_interface::nrg {

  /// The coefficients of one channel of a Wilson chain: the coupling θ and the on-site energies ξ and hoppings ζ
  struct channel_coefficients_t {
    double theta = 0;
    std::vector<double> xi = {}, zeta = {};
    bool operator==(channel_coefficients_t const &) const = default;
  };

  /// The coefficients of the Wilson chain of one twist, per coefficient channel (adapter A1: uncoupled channels)
  struct chain_coefficients_t {
    double z = 0, Lambda = 0, bandrescale = 1;
    std::vector<channel_coefficients_t> channels = {};
    bool operator==(chain_coefficients_t const &) const = default;
  };

  void h5_write(h5::group g, std::string const &name, chain_coefficients_t const &c);
  void h5_read(h5::group g, std::string const &name, chain_coefficients_t &c);

  /// Γ of the coefficient channels at the mesh points: gamma[k][c], diagonal in the channels
  struct channel_gamma_t {
    std::vector<double> omega;
    std::vector<std::vector<double>> gamma;
  };

  /// The keys of the [param] section of a parameter file: key -> value
  using param_section_t = std::map<std::string, std::string>;

  /// Nmax of a [param] section (resolve_chain_length: from Tmin and bandrescale, or Nmax as given)
  unsigned chain_length(param_section_t const &param);

  /// Static checks of a [param] section for the mixchain route (validate_preparation_parameters)
  void validate_chain_parameters(param_section_t const &param);

  /**
   * The chain coefficients of each twist: discretize_hybridization of Γ (one chain per z, nmax sites) and the
   * uncoupled-channel adapter (diagonal_coefficients). The [param] section gives mixchain's keys. With a non-empty
   * diagnostic_file the log goes there. Throws std::runtime_error with NRG Ljubljana's message.
   */
  std::vector<chain_coefficients_t> discretize(channel_gamma_t const &gamma, param_section_t const &param, std::vector<double> const &twists,
                                               unsigned nmax, std::string const &diagnostic_file = {});

  /// Initialize the calculation directory (its param, the template directory) from the chain coefficients; the ground
  /// energy of the seed
  double instantiate(std::string const &directory, std::string const &template_directory, unsigned nmax, chain_coefficients_t const &coefficients,
                     std::string const &diagnostic_file = {});

  /// Solve the initialized calculation directory, with scratch files in scratch_directory (absolute)
  void solve(std::string const &directory, std::string const &scratch_directory, std::string const &diagnostic_file = {});

} // namespace nrgljubljana_interface::nrg
