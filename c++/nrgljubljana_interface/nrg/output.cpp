#include "./output.hpp"

#include <nrgljubljana/solver_output.hpp>

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace nrgljubljana_interface::nrg {

  namespace {

    // NRG Ljubljana's readers name the line; the source is put in front
    template <typename F> auto with_source(std::string_view source, F &&f) {
      try {
        return f();
      } catch (std::runtime_error const &e) {
        if (source.empty()) throw;
        throw std::runtime_error(std::string(source) + ": " + e.what());
      }
    }

    // The spectrum list of a type letter of the manifest
    std::string spectrum_list(char type) {
      switch (type) {
        case 'd': return "specd";
        case 's': return "specs";
        case 't': return "spect";
        case 'o': return "specot";
        case 'q': return "specq";
        default: throw std::invalid_argument(std::string("no spectrum list of type '") + type + "'");
      }
    }

  } // namespace

  step_table_t parse_step_table(std::string_view text, std::string_view source) {
    auto const t = with_source(source, [&] { return NRG::Workflow::read_step_table(std::string(text)); });
    step_table_t result;
    result.columns = t.columns;
    result.values  = nda::array<double, 2>(long(t.rows.size()), long(t.columns.size()));
    for (long i = 0; i < long(t.rows.size()); ++i)
      for (long j = 0; j < long(t.columns.size()); ++j) result.values(i, j) = t.rows[i][j];
    return result;
  }

  std::map<std::string, double> parse_final_values(std::string_view text, std::string_view source) {
    return with_source(source, [&] { return NRG::Workflow::read_final_values(std::string(text)); });
  }

  energy_flow_t parse_energy_flow(std::string_view text, std::string_view source) {
    auto const flow = with_source(source, [&] { return NRG::Workflow::read_energy_flow(std::string(text)); });
    energy_flow_t result(flow.size());
    for (std::size_t n = 0; n < flow.size(); ++n)
      for (auto const &level : flow[n]) result[n].push_back({level.energy, level.subspaces, level.degeneracy});
    return result;
  }

  std::string spectrum_file(std::string const &op1, std::string const &op2, std::string const &suffix, char type) {
    auto const list = spectrum_list(type);
    if (op1.empty() or op2.empty()) throw std::invalid_argument("spectrum_file: an empty operator name");
    if (not suffix.empty() and suffix != "u" and suffix != "d")
      throw std::invalid_argument("spectrum_file: the spin suffix must be empty, u or d, got \"" + suffix + "\"");
    int const spin = suffix.empty() ? 0 : (suffix == "u" ? 1 : -1);
    return NRG::Workflow::spectrum_file_name(list, op1 + "-" + op2, "FDM", spin) + ".bin";
  }

  std::string read_text_file(std::string const &path) {
    std::ifstream f(path, std::ios::binary);
    if (not f) throw std::runtime_error(path + ": cannot be read");
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
  }

} // namespace nrgljubljana_interface::nrg
