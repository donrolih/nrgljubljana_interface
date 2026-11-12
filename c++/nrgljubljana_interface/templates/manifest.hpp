// The manifest of a template: what the interface needs to know about a template beyond its generated files.

#pragma once

#include <compare>
#include <map>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace nrgljubljana_interface::templates {

  /// Block names and sizes; the same type as triqs::gfs::gf_struct_t
  using block_structure_t = std::vector<std::pair<std::string, long>>;

  /// A block element "name[i,j]", zero-based
  struct block_element_t {
    std::string block;
    long i = 0, j = 0;
    auto operator<=>(block_element_t const &) const = default;
  };

  /// Parse "name[i,j]"; throws std::invalid_argument
  block_element_t parse_block_element(std::string_view text);

  /// "name[i,j]"
  std::string to_string(block_element_t const &e);

  /// One block element of a logical correlator and the NRG correlator that gives it
  struct correlator_entry_t {
    block_element_t element;
    std::string nrg_name;  // e.g. "A_d-A_d-u"
    std::string op1, op2;  // its operators, e.g. "A_d" and "A_d"
    std::string suffix;    // "", "u" or "d": the spin of a polarized spectrum
  };

  /// Names of the logical correlators with a fixed meaning: the correlators of the self-energy trick
  inline constexpr char const *self_energy_correlators[] = {"G", "F_l", "F_r", "I"};

  struct manifest_t {
    // [template]
    int format = 0;
    std::string model, symtype;
    long channels = 0;
    std::vector<std::string> ops, params;

    // [structure]
    block_structure_t gf, delta;
    std::map<std::string, long> degeneracy;

    // [bath]
    std::string layout = "uncoupled";
    std::vector<block_element_t> coefficient_channels;

    // [correlators]: logical name -> entries, ordered by element
    std::map<std::string, std::vector<correlator_entry_t>> correlators;

    // [hartree], [levels]: element -> operator, element -> model parameter
    std::map<block_element_t, std::string> hartree, levels;

    /// Whether the template gives F_l, F_r and I with G, and so the self-energy
    [[nodiscard]] bool has_self_energy() const;

    /// The block structure of a logical correlator: gf for G, F_l, F_r and I; otherwise the blocks its entries
    /// name, in alphabetical order, each sized by the largest index used in it
    [[nodiscard]] block_structure_t correlator_structure(std::string const &name) const;
  };

  /**
   * Parse manifest text and check it on its own (without the template files).
   *
   * Errors (std::invalid_argument, prefixed by source): the syntax errors of read_parameter_sections; an unknown
   * section or key; a missing required section or key; an unknown format; a malformed value; a block or list entry
   * given twice; a block element outside its structure; degeneracy not covering exactly the blocks of gf; an
   * unknown layout; a correlator whose operators are not in ops; F_l, F_r, I not given together with G for the
   * same elements, or without a [hartree] entry for each of them; a [hartree] entry for another element; a
   * [levels] entry that is off-diagonal or names a parameter not in params.
   */
  manifest_t parse_manifest(std::string_view text, std::string_view source = {});

  /// Read and parse <template_dir>/manifest
  manifest_t load_manifest(std::string const &template_dir);

  /// The outcome of checking a manifest against its template's files
  struct template_check_t {
    /// The NRG spectrum list of each NRG correlator in [correlators], by its letter: s (specs), t (spect),
    /// o (specot), d (specd), q (specq)
    std::map<std::string, char> correlator_type;
    /// Whether the symmetry type has a field (doublet spectra per spin)
    bool field = false;
    /// Checks that could not be run, each with the reason
    std::vector<std::string> pending;
  };

  /**
   * Check a manifest against the generated files of its template, through NRG Ljubljana's inspect_template_bundle.
   *
   * The template must pass the inspection; symtype, channels and the set of operators must match the manifest.
   * The operator types of a correlator must form a spectrum list as in NRG Ljubljana: s-s, p-p, g-g, s-g, g-s
   * (specs), t-t (spect), o-o (specot), d-d (specd), q-q (specq). A d-d correlator carries the spin suffix -u/-d
   * exactly when the symmetry type has a field; other correlators never do. A [hartree] operator must be of type s.
   * The number of coefficient channels must be channels or twice that.
   *
   * Throws std::invalid_argument naming the template, the file and the entry on a mismatch.
   */
  template_check_t check_against_template(manifest_t const &m, std::string const &template_dir);

} // namespace nrgljubljana_interface::templates
