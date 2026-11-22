// Parameter-file text: sectioned key=value lines, read with NRG Ljubljana's reader.

#pragma once

#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace nrgljubljana_interface::nrg {

  /// The keys of one section: key -> value
  using parameter_section_t = std::map<std::string, std::string>;

  /// All sections of a parameter file: section name -> keys
  using parameter_sections_t = std::map<std::string, parameter_section_t>;

  /**
   * Parse parameter-file text with NRG Ljubljana's reader (read_parameter_sections of the workflow library), so that
   * it reads as the solver reads its parameter files.
   *
   * Syntax: each line is trimmed of surrounding whitespace. Blank lines and lines starting with '#' are skipped
   * (comments are whole lines only; a '#' after a value is part of the value). A section header is a line that is
   * exactly "[name]" (no whitespace around or inside the brackets). "key=value" is split at the first '='; key and
   * value are trimmed. Keys may contain spaces and brackets, values may be empty.
   *
   * Errors (std::invalid_argument, with the line number, prefixed by source if given): a line that is none of the
   * above, a malformed section header or one with whitespace, a section given twice, a key before the first section,
   * an empty key, and a key given twice in one section.
   */
  parameter_sections_t read_parameter_sections(std::string_view text, std::string_view source = {});

  /// Split a list value at whitespace: "U1 eps1" -> {"U1", "eps1"}; an empty value gives an empty list
  std::vector<std::string> split_list(std::string_view value);

} // namespace nrgljubljana_interface::nrg
