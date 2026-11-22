#include "./parameter_file.hpp"

#include <nrgljubljana/solver_preparation.hpp>

#include <cstddef>
#include <stdexcept>

namespace nrgljubljana_interface::nrg {

  parameter_sections_t read_parameter_sections(std::string_view text, std::string_view source) {
    try {
      return NRG::Workflow::read_parameter_sections(std::string(text));
    } catch (std::runtime_error const &e) {
      throw std::invalid_argument(source.empty() ? std::string(e.what()) : std::string(source) + ": " + e.what());
    }
  }

  std::vector<std::string> split_list(std::string_view value) {
    constexpr std::string_view whitespace = " \t\r\n\f\v";
    std::vector<std::string> items;
    std::size_t pos = 0;
    while ((pos = value.find_first_not_of(whitespace, pos)) != std::string_view::npos) {
      auto const end = value.find_first_of(whitespace, pos);
      items.emplace_back(value.substr(pos, end == std::string_view::npos ? std::string_view::npos : end - pos));
      pos = end;
    }
    return items;
  }

} // namespace nrgljubljana_interface::nrg
