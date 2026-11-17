#include "./process.hpp"

#include "../nrg/broadening.hpp"
#include "../nrg/convert.hpp"

#include <numbers>
#include <stdexcept>

namespace nrgljubljana_interface {

  namespace {

    // Mean over the twists of a map of values, without the temperature column
    std::map<std::string, double> twist_average(std::vector<twist_result_t> const &twists, std::map<std::string, double> twist_result_t::*field) {
      std::map<std::string, double> sum;
      for (auto const &t : twists)
        for (auto const &[key, value] : t.*field)
          if (key != "T") sum[key] += value;
      for (auto &[key, value] : sum) value /= double(twists.size());
      return sum;
    }

    g_w_t zero_gf(triqs::mesh::refreq_log const &mesh, templates::block_structure_t const &s) {
      g_w_t g{mesh, s};
      for (auto &b : g.data()) b.data() = 0;
      return g;
    }

    long block_index(templates::block_structure_t const &s, std::string const &block) {
      for (long b = 0; b < long(s.size()); ++b)
        if (s[b].first == block) return b;
      throw std::logic_error("block " + block + " not in the structure");
    }

  } // namespace

  void process(solve_result_t &r, templates::manifest_t const &m, triqs::mesh::refreq_log const &mesh) {
    r.expv  = twist_average(r.twists, &twist_result_t::customfdm);
    r.tdfdm = twist_average(r.twists, &twist_result_t::tdfdm);
    r.processed.clear();
    if (not r.request.broaden or r.twists.empty()) return;

    auto const points = nrg::mesh_points(mesh);
    nrg::broadening_t const b{r.solve_params.alpha, r.solve_params.T, r.solve_params.gamma};

    for (auto const &name : requested_correlators(r.request, m)) {
      if (name == "Sigma" or name == "SigmaHartree") throw std::invalid_argument("correlator name " + name + " is reserved");
      auto const s = m.correlator_structure(name);
      auto g       = zero_gf(mesh, s);
      for (auto const &entry : m.correlators.at(name)) {
        std::vector<nrg::spectrum_bins_t const *> twists;
        for (auto const &t : r.twists) twists.push_back(&t.spectra.at(entry.nrg_name));
        auto const density = nrg::broaden(twists, points, b);

        auto const type      = r.correlator_type.at(entry.nrg_name);
        bool const fermionic = type == 'd' or type == 'q';
        std::vector<double> im(density.size());
        for (std::size_t k = 0; k < density.size(); ++k) im[k] = fermionic ? -std::numbers::pi * density[k] : density[k];
        auto const re = nrg::kramers_kronig(points, im);

        auto &gb = g[block_index(s, entry.element.block)];
        nrg::from_table({points, re}, gb, entry.element.i, entry.element.j, nrg::part_t::real);
        nrg::from_table({points, im}, gb, entry.element.i, entry.element.j, nrg::part_t::imag);
      }
      r.processed[name] = std::move(g);
    }

    if (not(r.processed.contains("G") and r.processed.contains("F_l") and r.processed.contains("F_r") and r.processed.contains("I"))) return;

    // Hartree term: twist average of the expectation value of the [hartree] operator, constant in frequency
    auto hartree = zero_gf(mesh, m.gf);
    for (auto const &[e, op] : m.hartree) {
      if (not op) continue;
      double value = 0;
      for (auto const &t : r.twists) {
        auto const it = t.customfdm.find(*op);
        if (it == t.customfdm.end()) throw std::runtime_error("expectation value " + *op + " of [hartree] not in customfdm");
        value += it->second;
      }
      value /= double(r.twists.size());
      auto &hb = hartree[block_index(m.gf, e.block)];
      for (auto w : mesh) hb[w](e.i, e.j) = value;
    }

    auto const &G = r.processed.at("G");
    auto const &F_l = r.processed.at("F_l");
    auto const &F_r = r.processed.at("F_r");
    auto const &I   = r.processed.at("I");
    g_w_t sigma                 = hartree + I - F_l / G * F_r;
    r.processed["Sigma"]        = std::move(sigma);
    r.processed["SigmaHartree"] = std::move(hartree);
  }

} // namespace nrgljubljana_interface
