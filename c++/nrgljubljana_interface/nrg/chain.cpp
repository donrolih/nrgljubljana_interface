#include "./chain.hpp"

#include <nrgljubljana/solver_chain.hpp>
#include <nrgljubljana/solver_preparation.hpp>
#include <nrgljubljana/solver_workflow.hpp>

#include <stdexcept>

namespace nrgljubljana_interface::nrg {

  namespace {

    chain_coefficients_t convert(NRG::Workflow::ChainCoefficients const &c) {
      chain_coefficients_t result{c.z, c.Lambda, c.bandrescale, {}};
      for (auto const &ch : c.channels) result.channels.push_back({ch.theta, ch.xi, ch.zeta});
      return result;
    }

    NRG::Workflow::ChainCoefficients convert(chain_coefficients_t const &c) {
      NRG::Workflow::ChainCoefficients result;
      result.z           = c.z;
      result.Lambda      = c.Lambda;
      result.bandrescale = c.bandrescale;
      for (auto const &ch : c.channels) result.channels.push_back({ch.theta, ch.xi, ch.zeta});
      return result;
    }

  } // namespace

  void h5_write(h5::group g, std::string const &name, chain_coefficients_t const &c) {
    auto grp = g.create_group(name);
    h5::write(grp, "z", c.z);
    h5::write(grp, "Lambda", c.Lambda);
    h5::write(grp, "bandrescale", c.bandrescale);
    auto channels = grp.create_group("channels");
    for (std::size_t k = 0; k < c.channels.size(); ++k) {
      auto ch = channels.create_group(std::to_string(k));
      h5::write(ch, "theta", c.channels[k].theta);
      h5::write(ch, "xi", c.channels[k].xi);
      h5::write(ch, "zeta", c.channels[k].zeta);
    }
  }

  void h5_read(h5::group g, std::string const &name, chain_coefficients_t &c) {
    auto grp = g.open_group(name);
    h5::read(grp, "z", c.z);
    h5::read(grp, "Lambda", c.Lambda);
    h5::read(grp, "bandrescale", c.bandrescale);
    auto channels = grp.open_group("channels");
    c.channels.clear();
    for (std::size_t k = 0; channels.has_key(std::to_string(k)); ++k) {
      auto ch   = channels.open_group(std::to_string(k));
      auto &out = c.channels.emplace_back();
      h5::read(ch, "theta", out.theta);
      h5::read(ch, "xi", out.xi);
      h5::read(ch, "zeta", out.zeta);
    }
  }

  unsigned chain_length(param_section_t const &param) { return NRG::Workflow::resolve_chain_length(param); }

  void validate_chain_parameters(param_section_t const &param) {
    NRG::Workflow::validate_preparation_parameters(param, NRG::Workflow::ChainSource::mixchain);
  }

  std::vector<chain_coefficients_t> discretize(channel_gamma_t const &gamma, param_section_t const &param, std::vector<double> const &twists,
                                               unsigned nmax, std::string const &diagnostic_file) {
    if (gamma.gamma.size() != gamma.omega.size()) throw std::invalid_argument("discretize: Γ and the mesh differ in size");
    NRG::Workflow::HybridizationFunction h;
    h.channels   = gamma.gamma.empty() ? 0 : unsigned(gamma.gamma.front().size());
    h.omega      = gamma.omega;
    unsigned n   = h.channels;
    h.real.assign(gamma.omega.size() * n * n, 0.0);
    for (std::size_t k = 0; k < gamma.omega.size(); ++k) {
      if (gamma.gamma[k].size() != n) throw std::invalid_argument("discretize: Γ has a different number of channels at mesh point " + std::to_string(k));
      for (unsigned c = 0; c < n; ++c) h.real[(k * n + c) * n + c] = gamma.gamma[k][c];
    }
    auto const chains = NRG::Workflow::discretize_hybridization(h, param, twists, nmax, diagnostic_file);
    std::vector<chain_coefficients_t> result;
    for (auto const &chain : chains) result.push_back(convert(NRG::Workflow::diagonal_coefficients(chain)));
    return result;
  }

  double instantiate(std::string const &directory, std::string const &template_directory, unsigned nmax, chain_coefficients_t const &coefficients,
                     std::string const &diagnostic_file) {
    NRG::Workflow::InitializerOptions options;
    options.calculation_directory = directory;
    options.parameter_file        = "param";
    options.template_directory    = template_directory;
    options.nmax                  = nmax;
    return NRG::Workflow::instantiate(options, convert(coefficients), diagnostic_file).ground_energy;
  }

  void solve(std::string const &directory, std::string const &scratch_directory, std::string const &diagnostic_file) {
    if (diagnostic_file.empty())
      NRG::Workflow::solve(directory, scratch_directory);
    else
      NRG::Workflow::solve(directory, scratch_directory, diagnostic_file);
  }

} // namespace nrgljubljana_interface::nrg
