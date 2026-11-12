#include "./parameter_file.hpp"

#include <cstddef>
#include <stdexcept>

namespace nrgljubljana_interface::nrg {

  namespace {

    constexpr std::string_view whitespace = " \t\r\n\f\v";

    std::string_view trim(std::string_view s) {
      auto const first = s.find_first_not_of(whitespace);
      if (first == std::string_view::npos) return {};
      auto const last = s.find_last_not_of(whitespace);
      return s.substr(first, last - first + 1);
    }

    [[noreturn]] void fail(std::string_view source, std::size_t line_number, std::string_view line, std::string const &what) {
      std::string msg;
      if (not source.empty()) msg += std::string(source) + ": ";
      msg += "line " + std::to_string(line_number) + ": " + what + ": \"" + std::string(line) + "\"";
      throw std::invalid_argument(msg);
    }

  } // namespace

  parameter_sections_t read_parameter_sections(std::string_view text, std::string_view source) {
    parameter_sections_t sections;
    parameter_section_t *current = nullptr;
    std::size_t line_number      = 0;

    while (not text.empty()) {
      auto const eol = text.find('\n');
      auto const raw = text.substr(0, eol);
      text           = eol == std::string_view::npos ? std::string_view{} : text.substr(eol + 1);
      ++line_number;

      auto const line = trim(raw);
      if (line.empty() or line.front() == '#') continue;

      if (line.front() == '[') {
        if (line.back() != ']') fail(source, line_number, line, "malformed section header");
        auto const name = std::string(trim(line.substr(1, line.size() - 2)));
        if (name.empty()) fail(source, line_number, line, "empty section name");
        auto [it, inserted] = sections.try_emplace(name);
        if (not inserted) fail(source, line_number, line, "section [" + name + "] given twice");
        current = &it->second;
        continue;
      }

      auto const pos_eq = line.find('=');
      if (pos_eq == std::string_view::npos) fail(source, line_number, line, "expected key=value, [section] or a # comment");
      if (current == nullptr) fail(source, line_number, line, "key before the first section");
      auto const key = std::string(trim(line.substr(0, pos_eq)));
      if (key.empty()) fail(source, line_number, line, "empty key");
      if (not current->try_emplace(key, trim(line.substr(pos_eq + 1))).second) fail(source, line_number, line, "key \"" + key + "\" given twice");
    }
    return sections;
  }

  std::vector<std::string> split_list(std::string_view value) {
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
