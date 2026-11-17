#include "./broadening.hpp"
#include "./convert.hpp"

#include <nrgljubljana/numerics.hpp>

#include <stdexcept>

namespace nrgljubljana_interface::nrg {

  std::vector<double> broaden(std::vector<spectrum_bins_t const *> const &twists, std::vector<double> const &mesh, broadening_t const &b) {
    if (twists.empty()) throw std::invalid_argument("broaden: no twists");
    std::vector<NRG::Numerics::SpectrumBins> bins;
    bins.reserve(twists.size());
    for (auto const *t : twists) bins.push_back(to_nrg(*t, part_t::real));

    NRG::Numerics::BroadeningOptions options{};
    options.alpha             = b.alpha;
    options.temperature       = b.temperature;
    options.omega0_ratio      = b.omega0_ratio;
    options.final_derfd_gamma = b.gamma;
    auto const t = NRG::Numerics::broaden(bins, mesh, options);
    if (t.x != mesh) throw std::runtime_error("broaden: result is not on the requested mesh");
    return t.y;
  }

  std::vector<double> kramers_kronig(std::vector<double> const &mesh, std::vector<double> const &values) {
    auto const t = NRG::Numerics::kramers_kronig(NRG::Numerics::Table{mesh, values});
    if (t.x != mesh) throw std::runtime_error("kramers_kronig: result is not on the input mesh");
    return t.y;
  }

} // namespace nrgljubljana_interface::nrg
