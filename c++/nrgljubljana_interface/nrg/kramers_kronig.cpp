#include "./kramers_kronig.hpp"

#include <nrgljubljana/numerics.hpp>

#include <stdexcept>

namespace nrgljubljana_interface::nrg {

  std::vector<double> kramers_kronig(std::vector<double> const &mesh, std::vector<double> const &values) {
    auto const t = NRG::Numerics::kramers_kronig(NRG::Numerics::Table{mesh, values});
    if (t.x != mesh) throw std::runtime_error("kramers_kronig: result is not on the input mesh");
    return t.y;
  }

} // namespace nrgljubljana_interface::nrg
