#include "./manifest.hpp"

#include "../nrg/parameter_file.hpp"

#include <nrgljubljana/solver_preparation.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace nrgljubljana_interface::templates {

  namespace {

    using nrg::split_list;

    constexpr std::string_view whitespace = " \t\r\n\f\v";

    [[noreturn]] void fail(std::string_view source, std::string const &where, std::string const &what) {
      std::string msg;
      if (not source.empty()) msg += std::string(source) + ": ";
      if (not where.empty()) msg += where + ": ";
      throw std::invalid_argument(msg + what);
    }

    std::string where(std::string const &section, std::string const &key = {}) { return "[" + section + "]" + (key.empty() ? "" : " " + key); }

    // A non-negative integer, the whole of text; -1 otherwise
    long parse_count(std::string_view text) {
      if (text.empty() or text.size() > 9 or text.find_first_not_of("0123456789") != std::string_view::npos) return -1;
      return std::stol(std::string(text));
    }

    bool valid_block_name(std::string_view name) { return not name.empty() and name.find_first_of("[],:= \t") == std::string_view::npos; }

    bool contains(std::vector<std::string> const &v, std::string const &x) { return std::find(v.begin(), v.end(), x) != v.end(); }

    // "name:size name:size ..."
    block_structure_t parse_structure(std::string const &value, std::string_view source, std::string const &at) {
      block_structure_t result;
      for (auto const &item : split_list(value)) {
        auto const colon = item.find(':');
        auto const name  = item.substr(0, colon);
        auto const size  = colon == std::string::npos ? -1 : parse_count(std::string_view(item).substr(colon + 1));
        if (not valid_block_name(name) or size < 1) fail(source, at, "expected name:size with size >= 1, got \"" + item + "\"");
        for (auto const &[n, s] : result)
          if (n == name) fail(source, at, "block \"" + name + "\" given twice");
        result.emplace_back(name, size);
      }
      return result;
    }

    std::vector<std::string> parse_unique_list(std::string const &value, std::string_view source, std::string const &at) {
      auto result = split_list(value);
      for (std::size_t k = 0; k < result.size(); ++k)
        if (std::find(result.begin(), result.begin() + k, result[k]) != result.begin() + k) fail(source, at, "\"" + result[k] + "\" given twice");
      return result;
    }

    block_element_t parse_element(std::string_view text, std::string_view source, std::string const &at) {
      try {
        return parse_block_element(text);
      } catch (std::invalid_argument const &e) { fail(source, at, e.what()); }
    }

    void check_in_structure(block_element_t const &e, block_structure_t const &s, std::string const &s_name, std::string_view source,
                            std::string const &at) {
      for (auto const &[name, size] : s)
        if (name == e.block) {
          if (e.i >= size or e.j >= size) fail(source, at, to_string(e) + " is outside block \"" + name + "\" of size " + std::to_string(size));
          return;
        }
      fail(source, at, "no block \"" + e.block + "\" in " + s_name);
    }

    // Split an NRG correlator "op1-op2" or "op1-op2-u|d" into its operators; operator names may contain '-'
    void split_correlator(correlator_entry_t &entry, std::vector<std::string> const &ops, std::string_view source, std::string const &at) {
      auto const &name = entry.nrg_name;
      std::vector<std::array<std::string, 3>> candidates;
      for (auto pos = name.find('-'); pos != std::string::npos; pos = name.find('-', pos + 1)) {
        auto const left = name.substr(0, pos), rest = name.substr(pos + 1);
        if (not contains(ops, left)) continue;
        if (contains(ops, rest)) candidates.push_back({left, rest, ""});
        if (rest.size() > 2 and rest[rest.size() - 2] == '-' and (rest.back() == 'u' or rest.back() == 'd')
            and contains(ops, rest.substr(0, rest.size() - 2)))
          candidates.push_back({left, rest.substr(0, rest.size() - 2), rest.substr(rest.size() - 1)});
      }
      if (candidates.empty()) fail(source, at, "\"" + name + "\" is not <op>-<op> or <op>-<op>-u|d with operators in ops");
      if (candidates.size() > 1) fail(source, at, "\"" + name + "\" splits into operators of ops in more than one way");
      entry.op1    = candidates[0][0];
      entry.op2    = candidates[0][1];
      entry.suffix = candidates[0][2];
    }

    std::set<block_element_t> elements_of(std::vector<correlator_entry_t> const &entries) {
      std::set<block_element_t> result;
      for (auto const &e : entries) result.insert(e.element);
      return result;
    }

    // The spectrum list (s, t, o, d, q) of a correlator of operator types t1 and t2, or 0 if there is none
    char spectrum_list(char t1, char t2) {
      auto const singlet = [](char t) { return t == 's' or t == 'p' or t == 'g'; };
      if (singlet(t1) and singlet(t2)) return (t1 == t2 or (t1 != 'p' and t2 != 'p')) ? 's' : 0;
      if (t1 != t2) return 0;
      return std::string_view("tdoq").find(t1) == std::string_view::npos ? 0 : t1;
    }

    std::string read_file(std::filesystem::path const &path) {
      std::ifstream f(path, std::ios::binary);
      if (not f) throw std::invalid_argument(path.string() + ": cannot be read");
      std::stringstream ss;
      ss << f.rdbuf();
      return ss.str();
    }

  } // namespace

  // ---------------------------------------------------------------------------------------------------------------

  block_element_t parse_block_element(std::string_view text) {
    auto const bad = [&]() -> block_element_t { throw std::invalid_argument("expected a block element name[i,j], got \"" + std::string(text) + "\""); };
    auto const open = text.find('[');
    if (open == std::string_view::npos or text.back() != ']') return bad();
    auto const name  = text.substr(0, open);
    auto const inner = text.substr(open + 1, text.size() - open - 2);
    auto const comma = inner.find(',');
    if (not valid_block_name(name) or comma == std::string_view::npos) return bad();
    auto const i = parse_count(inner.substr(0, comma)), j = parse_count(inner.substr(comma + 1));
    if (i < 0 or j < 0) return bad();
    return {std::string(name), i, j};
  }

  std::string to_string(block_element_t const &e) { return e.block + "[" + std::to_string(e.i) + "," + std::to_string(e.j) + "]"; }

  bool has_gf_structure(std::string_view name) {
    return std::any_of(std::begin(gf_correlators), std::end(gf_correlators), [&](char const *s) { return name == s; });
  }

  bool manifest_t::has_self_energy() const { return correlators.contains("F_l"); }

  block_structure_t manifest_t::correlator_structure(std::string const &name) const {
    if (has_gf_structure(name)) return gf;
    std::map<std::string, long> sizes;
    if (auto it = correlators.find(name); it != correlators.end())
      for (auto const &entry : it->second) {
        auto &size = sizes[entry.element.block];
        size       = std::max({size, entry.element.i + 1, entry.element.j + 1});
      }
    return {sizes.begin(), sizes.end()};
  }

  // ---------------------------------------------------------------------------------------------------------------

  manifest_t parse_manifest(std::string_view text, std::string_view source) {
    auto sections = nrg::read_parameter_sections(text, source);

    // Known sections and keys; an empty key list means any key
    std::map<std::string, std::vector<std::string>> const known = {
       {"template", {"format", "model", "symtype", "channels", "ops", "params"}},
       {"structure", {"gf", "delta", "degeneracy"}},
       {"bath", {"layout", "coefficient_channels"}},
       {"correlators", {}},
       {"hartree", {}},
       {"levels", {}},
    };
    for (auto const &[section, keys] : sections) {
      auto const it = known.find(section);
      if (it == known.end()) fail(source, "", "unknown section [" + section + "]");
      if (it->second.empty()) continue;
      for (auto const &[key, value] : keys)
        if (not contains(it->second, key)) fail(source, where(section), "unknown key \"" + key + "\"");
    }

    auto const required = [&](std::string const &section, std::string const &key) -> std::string const & {
      auto const s = sections.find(section);
      if (s == sections.end()) fail(source, "", "missing section [" + section + "]");
      auto const k = s->second.find(key);
      if (k == s->second.end()) fail(source, where(section), "missing key \"" + key + "\"");
      return k->second;
    };
    auto const section = [&](std::string const &name) -> nrg::parameter_section_t const & {
      static nrg::parameter_section_t const empty;
      auto const s = sections.find(name);
      return s == sections.end() ? empty : s->second;
    };

    manifest_t m;

    // [template]
    if (required("template", "format") != "1") fail(source, where("template", "format"), "unknown format \"" + required("template", "format") + "\"");
    m.format  = 1;
    m.model   = required("template", "model");
    m.symtype = required("template", "symtype");
    if (m.model.empty()) fail(source, where("template", "model"), "empty");
    if (m.symtype.empty()) fail(source, where("template", "symtype"), "empty");
    m.channels = parse_count(required("template", "channels"));
    if (m.channels < 1) fail(source, where("template", "channels"), "expected a positive integer");
    m.ops    = parse_unique_list(required("template", "ops"), source, where("template", "ops"));
    m.params = parse_unique_list(required("template", "params"), source, where("template", "params"));

    // [structure]
    m.gf = parse_structure(required("structure", "gf"), source, where("structure", "gf"));
    if (section("structure").contains("delta"))
      m.delta = parse_structure(section("structure").at("delta"), source, where("structure", "delta"));
    else
      m.delta = m.gf;
    if (m.delta.empty()) fail(source, where("structure", "delta"), "no blocks of Δ (gf is empty, so delta must be given)");
    for (auto const &[name, d] : parse_structure(required("structure", "degeneracy"), source, where("structure", "degeneracy"))) {
      if (std::none_of(m.gf.begin(), m.gf.end(), [&](auto const &b) { return b.first == name; }))
        fail(source, where("structure", "degeneracy"), "no block \"" + name + "\" in gf");
      m.degeneracy[name] = d;
    }
    for (auto const &[name, size] : m.gf)
      if (not m.degeneracy.contains(name)) fail(source, where("structure", "degeneracy"), "missing for block \"" + name + "\"");

    // [bath]
    if (section("bath").contains("layout")) m.layout = section("bath").at("layout");
    if (m.layout != "uncoupled" and m.layout != "pol2x2" and m.layout != "nambu")
      fail(source, where("bath", "layout"), "unknown layout \"" + m.layout + "\" (uncoupled, pol2x2 or nambu)");
    for (auto const &item : parse_unique_list(required("bath", "coefficient_channels"), source, where("bath", "coefficient_channels"))) {
      auto const e = parse_element(item, source, where("bath", "coefficient_channels"));
      check_in_structure(e, m.delta, "delta", source, where("bath", "coefficient_channels"));
      m.coefficient_channels.push_back(e);
    }
    if (m.coefficient_channels.empty()) fail(source, where("bath", "coefficient_channels"), "empty");

    // [correlators]
    for (auto const &[key, value] : section("correlators")) {
      auto const at = where("correlators", key);
      auto const f  = split_list(key);
      if (f.size() != 2) fail(source, at, "expected \"<logical name> <block element>\"");
      correlator_entry_t entry{parse_element(f[1], source, at), value, {}, {}, {}};
      if (value.empty() or value.find_first_of(whitespace) != std::string::npos) fail(source, at, "expected one NRG correlator, got \"" + value + "\"");
      split_correlator(entry, m.ops, source, at);
      auto const &name = f[0];
      if (has_gf_structure(name)) check_in_structure(entry.element, m.gf, "gf", source, at);
      auto &entries = m.correlators[name];
      for (auto const &other : entries)
        if (other.element == entry.element) fail(source, at, to_string(entry.element) + " of " + name + " given twice");
      entries.push_back(entry);
    }
    for (auto &[name, entries] : m.correlators)
      std::sort(entries.begin(), entries.end(), [](auto const &a, auto const &b) { return a.element < b.element; });

    // The self-energy trick: F_l, F_r and I come with G, for the same elements, each with a Hartree term
    bool const any_f = m.correlators.contains("F_l") or m.correlators.contains("F_r") or m.correlators.contains("I");
    std::set<block_element_t> sigma_elements;
    if (any_f) {
      if (not m.correlators.contains("G")) fail(source, where("correlators"), "F_l, F_r and I are given without G");
      sigma_elements = elements_of(m.correlators.at("G"));
      for (auto const *name : {"F_l", "F_r", "I"}) {
        if (not m.correlators.contains(name)) fail(source, where("correlators"), std::string(name) + " is missing; F_l, F_r and I are given together");
        if (elements_of(m.correlators.at(name)) != sigma_elements)
          fail(source, where("correlators"), std::string(name) + " is not given for the same elements as G");
      }
    }

    // [hartree]
    for (auto const &[key, value] : section("hartree")) {
      auto const at = where("hartree", key);
      auto const e  = parse_element(key, source, at);
      if (not sigma_elements.contains(e)) fail(source, at, "no self-energy for " + key + " (G, F_l, F_r and I give it)");
      if (value == "0") {
        m.hartree[e] = std::nullopt;
        continue;
      }
      if (not contains(m.ops, value)) fail(source, at, "operator \"" + value + "\" is not in ops (or 0 for no Hartree term)");
      m.hartree[e] = value;
    }
    for (auto const &e : sigma_elements)
      if (not m.hartree.contains(e)) fail(source, where("hartree"), "missing for " + to_string(e));

    // [levels]
    for (auto const &[key, value] : section("levels")) {
      auto const at = where("levels", key);
      auto const e  = parse_element(key, source, at);
      check_in_structure(e, m.gf, "gf", source, at);
      if (e.i != e.j) fail(source, at, "an impurity level belongs to a diagonal element");
      if (not contains(m.params, value)) fail(source, at, "parameter \"" + value + "\" is not in params");
      m.levels[e] = value;
    }

    return m;
  }

  manifest_t load_manifest(std::string const &template_dir) {
    auto const path = std::filesystem::path(template_dir) / "manifest";
    return parse_manifest(read_file(path), path.string());
  }

  // ---------------------------------------------------------------------------------------------------------------

  template_check_t check_against_template(manifest_t const &m, std::string const &template_dir) {
    auto const dir      = std::filesystem::path(template_dir);
    auto const manifest = (dir / "manifest").string();
    auto const data_in  = (dir / "data.in").string();

    std::map<std::string, std::string> files;
    for (auto const &f : std::filesystem::recursive_directory_iterator(dir))
      if (f.is_regular_file()) files.emplace(std::filesystem::relative(f.path(), dir).generic_string(), read_file(f.path()));
    NRG::Workflow::TemplateInfo info;
    try {
      info = NRG::Workflow::inspect_template_bundle(files);
    } catch (std::exception const &e) { fail(dir.string(), "inspect_template_bundle", e.what()); }

    template_check_t result;
    result.field        = info.field;
    auto const mismatch = [&](std::string const &at, std::string const &what) { fail(manifest, at, what + " (" + data_in + ")"); };

    if (m.symtype != info.symmetry) mismatch(where("template", "symtype"), "\"" + m.symtype + "\", but data.in states \"" + info.symmetry + "\"");
    if (m.channels != static_cast<long>(info.channels))
      mismatch(where("template", "channels"), std::to_string(m.channels) + ", but data.in states " + std::to_string(info.channels));

    // Operators: names and types
    std::map<std::string, char> type_of;
    for (auto const &op : info.operators) {
      auto const [it, inserted] = type_of.emplace(op.name, op.type);
      if (not inserted and it->second != op.type) fail(data_in, "", "operator \"" + op.name + "\" given with two types");
    }
    for (auto const &op : m.ops)
      if (not type_of.contains(op)) mismatch(where("template", "ops"), "operator \"" + op + "\" is not in data.in");
    for (auto const &[name, type] : type_of)
      if (not contains(m.ops, name)) mismatch(where("template", "ops"), "operator \"" + name + "\" of data.in is missing");

    for (auto const &[name, entries] : m.correlators)
      for (auto const &entry : entries) {
        auto const at   = where("correlators", name + " " + to_string(entry.element));
        auto const t1   = type_of.at(entry.op1), t2 = type_of.at(entry.op2);
        auto const list = spectrum_list(t1, t2);
        if (list == 0)
          mismatch(at, "operators \"" + entry.op1 + "\" (" + t1 + ") and \"" + entry.op2 + "\" (" + t2 + ") form no NRG spectrum");
        if (list == 'd' and info.field and entry.suffix.empty())
          mismatch(at, "the symmetry type " + info.symmetry + " has a field: a d spectrum needs the spin suffix -u or -d");
        if (not entry.suffix.empty() and not(list == 'd' and info.field))
          mismatch(at, "a spin suffix -u/-d belongs only to d spectra of a symmetry type with a field");
        result.correlator_type[entry.nrg_name] = list;
      }

    for (auto const &[e, op] : m.hartree)
      if (op and type_of.at(*op) != 's') mismatch(where("hartree", to_string(e)), "operator \"" + *op + "\" is not of type s");

    auto const n_coef = static_cast<long>(m.coefficient_channels.size());
    if (m.layout == "uncoupled" and n_coef != m.channels and n_coef != 2 * m.channels)
      mismatch(where("bath", "coefficient_channels"),
               std::to_string(n_coef) + " entries, but data.in has " + std::to_string(m.channels) + " channels (expected that or twice that)");
    result.pending.emplace_back("[bath] coefficient_channels: whether the run may be spin-polarized is not checked");
    if (m.layout != "uncoupled") result.pending.emplace_back("[bath] layout=" + m.layout + ": not checked against the template");
    result.pending.emplace_back("[template] params: not compared with the parameters used in the template's expressions");

    return result;
  }

} // namespace nrgljubljana_interface::templates
