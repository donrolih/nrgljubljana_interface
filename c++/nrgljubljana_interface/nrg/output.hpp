// Output files of an NRG Ljubljana run: tables of expectation values and thermodynamics, the energy flow, and the
// names of the spectrum files. The parsing is NRG Ljubljana's (solver_output.hpp); these give the interface's types.

#pragma once

#include <nda/nda.hpp>

#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace nrgljubljana_interface::nrg {

  /// A table with one row per NRG step (custom, td): column names and values
  struct step_table_t {

    /// Column names
    std::vector<std::string> columns;

    /// Values, indexed (step, column)
    nda::array<double, 2> values;
  };

  /// A group of degenerate levels of one NRG step, as in annotated.dat
  struct energy_level_t {
    /// Energy of the group
    double energy = 0;

    /// Quantum numbers of each invariant subspace in the group
    std::vector<std::vector<int>> subspaces;

    /// Total number of states, including multiplicity
    long degeneracy = 0;

    bool operator==(energy_level_t const &) const = default;
  };

  /// Lowest levels of each NRG step
  using energy_flow_t = std::vector<std::vector<energy_level_t>>;

  /**
   * Parse a table of custom, td, customfdm or tdfdm.
   *
   * The column names are those of the last '#' line before the data ("# T name ..."; custom has a line of column
   * numbers before it). Each further line is one row with a value per column.
   * Throws std::runtime_error, with the source and line, on malformed text.
   */
  step_table_t parse_step_table(std::string_view text, std::string_view source = {});

  /// Parse customfdm or tdfdm: a table with exactly one row, as column name -> value; errors as parse_step_table
  std::map<std::string, double> parse_final_values(std::string_view text, std::string_view source = {});

  /**
   * Parse annotated.dat with grouped levels: per line "energy (qn ...) (qn ...) ... [degeneracy]", one group of
   * degenerate levels; an empty line ends each NRG step. A file written with dumpgroups=false is refused.
   * Throws std::runtime_error, with the source and line, on malformed text.
   */
  energy_flow_t parse_energy_flow(std::string_view text, std::string_view source = {});

  /**
   * The file NRG Ljubljana writes the FDM spectral bins of an NRG correlator to, given its spectrum list (s, t, o, d,
   * q): "<prefix>_FDM_dens_<op1>-<op2>[-u|-d].bin" with the prefix corr, spin, orbspin, spec or specq. The operators
   * and the spin suffix ("", "u" or "d") are those the manifest splits the correlator into (operator names may
   * contain '-'). Throws std::invalid_argument for another type or suffix, or an empty operator name.
   */
  std::string spectrum_file(std::string const &op1, std::string const &op2, std::string const &suffix, char type);

  /// The contents of a file; throws std::runtime_error if it cannot be read
  std::string read_text_file(std::string const &path);

} // namespace nrgljubljana_interface::nrg
