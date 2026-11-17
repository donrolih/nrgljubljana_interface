// Template manifests (templates/manifest.hpp)

#include <gtest/gtest.h>

#include <nrgljubljana_interface/templates/manifest.hpp>

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace nrgljubljana_interface;
using namespace nrgljubljana_interface::templates;

namespace {

  std::string const siam_qs = R"(# SIAM, QS
[template]
format=1
model=SIAM
symtype=QS
channels=1
ops=n_d n_d^2 hop0 A_d self_d sigma_d SigmaHartree
params=U1 eps1

[structure]
gf=imp:1
degeneracy=imp:2

[bath]
coefficient_channels=imp[0,0]

[correlators]
G imp[0,0]=A_d-A_d
F_l imp[0,0]=self_d-A_d
F_r imp[0,0]=A_d-self_d
I imp[0,0]=self_d-self_d
chi_NN imp[0,0]=n_d-n_d
chi_SS imp[0,0]=sigma_d-sigma_d

[hartree]
imp[0,0]=SigmaHartree

[levels]
imp[0,0]=eps1
)";

  // text with the first occurrence of from replaced by to
  std::string replaced(std::string text, std::string const &from, std::string const &to) {
    auto const pos = text.find(from);
    if (pos == std::string::npos) throw std::logic_error("test text does not contain \"" + from + "\"");
    return text.replace(pos, from.size(), to);
  }

  // The error message of parse_manifest(text, "manifest"), or "" if it parses
  std::string error_of(std::string const &text) {
    try {
      parse_manifest(text, "manifest");
    } catch (std::invalid_argument const &e) { return e.what(); }
    return {};
  }

  std::string template_dir(std::string const &name) {
    auto const *root = std::getenv("NRGIF_TEMPLATE_DIR");
    if (root == nullptr) throw std::runtime_error("NRGIF_TEMPLATE_DIR is not set");
    return std::string(root) + "/" + name;
  }

  // The error message of check_against_template for a template directory, or "" if it passes
  std::string check_error_of(std::string const &text, std::string const &dir) {
    try {
      check_against_template(parse_manifest(text, "manifest"), dir);
    } catch (std::invalid_argument const &e) { return e.what(); }
    return {};
  }

  std::string check_error_of(std::string const &text) { return check_error_of(text, template_dir("SIAM/QS")); }

  // A copy of a shipped template in the working directory, with the first occurrence of from in data.in replaced by to
  std::string modified_template(std::string const &name, std::string const &copy, std::string const &from, std::string const &to) {
    namespace fs = std::filesystem;
    fs::remove_all(copy);
    fs::copy(template_dir(name), copy, fs::copy_options::recursive);
    std::stringstream ss;
    ss << std::ifstream(copy + "/data.in").rdbuf();
    std::ofstream(copy + "/data.in") << replaced(ss.str(), from, to);
    return copy;
  }

  std::string const siam_qsz = R"([template]
format=1
model=SIAM
symtype=QSZ
channels=1
ops=n_d n_d^2 hop0 A_d self_d SZd n_d_u n_d_d SigmaHartree-u SigmaHartree-d
params=U1 eps1 B1
[structure]
gf=up:1 dn:1
degeneracy=up:1 dn:1
[bath]
coefficient_channels=up[0,0] dn[0,0]
[correlators]
G up[0,0]=A_d-A_d-u
G dn[0,0]=A_d-A_d-d
chi_NN imp[0,0]=n_d-n_d
chi_SzSz imp[0,0]=SZd-SZd
)";

  bool contains(std::string const &s, std::string const &part) { return s.find(part) != std::string::npos; }

} // namespace

// ---- Block elements ----

TEST(Manifest, BlockElement) {
  auto const e = parse_block_element("up[1,0]");
  EXPECT_EQ(e.block, "up");
  EXPECT_EQ(e.i, 1);
  EXPECT_EQ(e.j, 0);
  EXPECT_EQ(to_string(e), "up[1,0]");
  for (auto const *bad : {"up", "up[1]", "up[1,0", "[0,0]", "up[-1,0]", "up[a,0]", "up[0,0]x", "up[ 0,0]", "u p[0,0]"})
    EXPECT_THROW(parse_block_element(bad), std::invalid_argument) << bad;
}

// ---- Valid manifests ----

TEST(Manifest, SiamQs) {
  auto const m = parse_manifest(siam_qs);
  EXPECT_EQ(m.format, 1);
  EXPECT_EQ(m.model, "SIAM");
  EXPECT_EQ(m.symtype, "QS");
  EXPECT_EQ(m.channels, 1);
  EXPECT_EQ(m.ops.size(), 7u);
  EXPECT_EQ(m.params, (std::vector<std::string>{"U1", "eps1"}));
  EXPECT_EQ(m.gf, (block_structure_t{{"imp", 1}}));
  EXPECT_EQ(m.delta, m.gf); // default
  EXPECT_EQ(m.degeneracy.at("imp"), 2);
  EXPECT_EQ(m.layout, "uncoupled"); // default
  ASSERT_EQ(m.coefficient_channels.size(), 1u);
  EXPECT_EQ(to_string(m.coefficient_channels[0]), "imp[0,0]");
  EXPECT_EQ(m.correlators.size(), 6u);
  auto const &f_l = m.correlators.at("F_l").at(0);
  EXPECT_EQ(f_l.nrg_name, "self_d-A_d");
  EXPECT_EQ(f_l.op1, "self_d");
  EXPECT_EQ(f_l.op2, "A_d");
  EXPECT_EQ(f_l.suffix, "");
  EXPECT_TRUE(m.has_self_energy());
  EXPECT_EQ(m.hartree.at(parse_block_element("imp[0,0]")), "SigmaHartree");
  EXPECT_EQ(m.levels.at(parse_block_element("imp[0,0]")), "eps1");
  EXPECT_EQ(m.correlator_structure("G"), m.gf);
  EXPECT_EQ(m.correlator_structure("chi_NN"), (block_structure_t{{"imp", 1}}));
}

TEST(Manifest, PolarizedAndDashedOperatorNames) {
  auto const m = parse_manifest(R"([template]
format=1
model=SIAM
symtype=QSZ
channels=1
ops=n_d A_d self_d SigmaHartree-u SigmaHartree-d
params=U1 eps1 B1
[structure]
gf=up:1 dn:1
degeneracy=up:1 dn:1
[bath]
coefficient_channels=up[0,0] dn[0,0]
[correlators]
G up[0,0]=A_d-A_d-u
G dn[0,0]=A_d-A_d-d
F_l up[0,0]=self_d-A_d-u
F_l dn[0,0]=self_d-A_d-d
F_r up[0,0]=A_d-self_d-u
F_r dn[0,0]=A_d-self_d-d
I up[0,0]=self_d-self_d-u
I dn[0,0]=self_d-self_d-d
chi_NN imp[0,0]=n_d-n_d
[hartree]
up[0,0]=SigmaHartree-u
dn[0,0]=SigmaHartree-d
[levels]
up[0,0]=eps1
dn[0,0]=eps1
)");
  auto const &g = m.correlators.at("G");
  ASSERT_EQ(g.size(), 2u);
  EXPECT_EQ(to_string(g[0].element), "dn[0,0]"); // ordered by element
  EXPECT_EQ(g[0].op1, "A_d");
  EXPECT_EQ(g[0].op2, "A_d");
  EXPECT_EQ(g[0].suffix, "d");
  EXPECT_EQ(g[1].suffix, "u");
  EXPECT_EQ(m.gf, (block_structure_t{{"up", 1}, {"dn", 1}})); // in the given order
  EXPECT_EQ(m.correlator_structure("chi_NN"), (block_structure_t{{"imp", 1}}));
  EXPECT_EQ(m.hartree.at(parse_block_element("dn[0,0]")), "SigmaHartree-d");
}

TEST(Manifest, CorrelatorStructureFromEntries) {
  auto const text = replaced(siam_qs, "chi_NN imp[0,0]=n_d-n_d", "chi_NN imp[0,0]=n_d-n_d\nchi_NN imp[0,2]=n_d-n_d\nchi_NN x[0,0]=n_d-n_d");
  auto const m = parse_manifest(text);
  EXPECT_EQ(m.correlator_structure("chi_NN"), (block_structure_t{{"imp", 3}, {"x", 1}}));
  EXPECT_EQ(m.correlator_structure("none"), block_structure_t{});
}

TEST(Manifest, GOnlyAndNoGf) {
  // G without F_l, F_r, I: no self-energy, no [hartree]
  std::string text = siam_qs;
  for (auto const *line : {"F_l imp[0,0]=self_d-A_d\n", "F_r imp[0,0]=A_d-self_d\n", "I imp[0,0]=self_d-self_d\n", "[hartree]\nimp[0,0]=SigmaHartree\n"})
    text = replaced(text, line, "");
  EXPECT_FALSE(parse_manifest(text).has_self_energy());

  // A bath-only template: no G blocks, Δ given
  auto const clean = parse_manifest(R"([template]
format=1
model=CLEAN
symtype=QS
channels=1
ops=n_f
params=
[structure]
gf=
delta=imp:1
degeneracy=
[bath]
coefficient_channels=imp[0,0]
)");
  EXPECT_TRUE(clean.gf.empty());
  EXPECT_EQ(clean.delta, (block_structure_t{{"imp", 1}}));
  EXPECT_TRUE(clean.params.empty());
}

// ---- Invalid manifests ----

TEST(Manifest, UnknownAndMissing) {
  EXPECT_TRUE(contains(error_of(siam_qs + "[level]\nimp[0,0]=eps1\n"), "unknown section [level]"));
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "model=SIAM", "modl=SIAM")), "[template]: unknown key \"modl\""));
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "[bath]\ncoefficient_channels=imp[0,0]\n", "")), "missing section [bath]"));
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "params=U1 eps1\n", "")), "[template]: missing key \"params\""));
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "degeneracy=imp:2\n", "")), "[structure]: missing key \"degeneracy\""));
  // syntax errors of the parameter-file parser come through with the source
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "model=SIAM", "model SIAM")), "manifest: line 4"));
}

TEST(Manifest, TemplateValues) {
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "format=1", "format=2")), "[template] format: unknown format \"2\""));
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "channels=1", "channels=0")), "[template] channels: expected a positive integer"));
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "channels=1", "channels=x")), "[template] channels"));
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "hop0 A_d", "hop0 hop0 A_d")), "[template] ops: \"hop0\" given twice"));
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "symtype=QS", "symtype=")), "[template] symtype: empty"));
}

TEST(Manifest, StructureValues) {
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "gf=imp:1", "gf=imp")), "[structure] gf: expected name:size"));
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "gf=imp:1", "gf=imp:0")), "[structure] gf: expected name:size"));
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "gf=imp:1", "gf=imp:1 imp:2")), "block \"imp\" given twice"));
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "degeneracy=imp:2", "degeneracy=imp:2 x:1")), "[structure] degeneracy: no block \"x\" in gf"));
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "gf=imp:1", "gf=imp:1 x:1")), "[structure] degeneracy: missing for block \"x\""));
  EXPECT_TRUE(contains(error_of(replaced(replaced(siam_qs, "gf=imp:1", "gf="), "degeneracy=imp:2", "degeneracy=")), "[structure] delta: no blocks"));
}

TEST(Manifest, BathValues) {
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "[bath]\n", "[bath]\nlayout=coupled\n")), "[bath] layout: unknown layout \"coupled\""));
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "coefficient_channels=imp[0,0]", "coefficient_channels=imp[1,1]")),
                       "[bath] coefficient_channels: imp[1,1] is outside block \"imp\" of size 1"));
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "coefficient_channels=imp[0,0]", "coefficient_channels=x[0,0]")), "no block \"x\" in delta"));
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "coefficient_channels=imp[0,0]", "coefficient_channels=imp[0,0] imp[0,0]")), "given twice"));
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "coefficient_channels=imp[0,0]", "coefficient_channels=")), "[bath] coefficient_channels: empty"));
}

TEST(Manifest, CorrelatorEntries) {
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "chi_NN imp[0,0]", "chi_NN")), "[correlators] chi_NN: expected \"<correlator> <block element>\""));
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "chi_NN imp[0,0]", "chi_NN imp(0,0)")), "expected a block element"));
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "=n_d-n_d", "=n_d-n_x")), "\"n_d-n_x\" is not <op>-<op>"));
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "=n_d-n_d", "=n_d n_d")), "expected one NRG correlator"));
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "G imp[0,0]", "G imp[0,1]")), "[correlators] G imp[0,1]: imp[0,1] is outside block \"imp\""));
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "chi_NN imp[0,0]=n_d-n_d", "chi_NN imp[0,0]=n_d-n_d\nchi_NN  imp[0,0]=n_d-n_d")),
                       "imp[0,0] of chi_NN given twice"));
}

TEST(Manifest, AmbiguousCorrelator) {
  // "a-b-u" splits as a + b-u and as a + b with suffix u
  auto const text = replaced(replaced(siam_qs, "sigma_d SigmaHartree", "sigma_d SigmaHartree a b b-u"), "=n_d-n_d", "=a-b-u");
  EXPECT_TRUE(contains(error_of(text), "splits into operators of ops in more than one way"));
}

TEST(Manifest, SelfEnergyRule) {
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "I imp[0,0]=self_d-self_d\n", "")), "I is missing"));
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "G imp[0,0]=A_d-A_d\n", "")), "F_l, F_r and I are given without G"));
  auto const two = replaced(replaced(siam_qs, "gf=imp:1", "gf=imp:2"), "G imp[0,0]=A_d-A_d", "G imp[0,0]=A_d-A_d\nG imp[1,1]=A_d-A_d");
  EXPECT_TRUE(contains(error_of(two), "F_l is not given for the same elements as G"));
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "[hartree]\nimp[0,0]=SigmaHartree\n", "")), "[hartree]: missing for imp[0,0]"));
}

TEST(Manifest, NoHartreeTerm) {
  auto const m = parse_manifest(replaced(siam_qs, "imp[0,0]=SigmaHartree", "imp[0,0]=0"));
  EXPECT_TRUE(m.has_self_energy());
  EXPECT_EQ(m.hartree.at(parse_block_element("imp[0,0]")), std::nullopt);
  EXPECT_EQ(check_error_of(replaced(siam_qs, "imp[0,0]=SigmaHartree", "imp[0,0]=0")), "");
}

TEST(Manifest, HartreeAndLevels) {
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "imp[0,0]=SigmaHartree", "imp[0,0]=SigmaH")), "[hartree] imp[0,0]: operator \"SigmaH\" is not in ops"));
  auto const g_only = replaced(replaced(replaced(siam_qs, "F_l imp[0,0]=self_d-A_d\n", ""), "F_r imp[0,0]=A_d-self_d\n", ""), "I imp[0,0]=self_d-self_d\n", "");
  EXPECT_TRUE(contains(error_of(g_only), "[hartree] imp[0,0]: no self-energy for imp[0,0]"));
  auto const two = replaced(siam_qs, "gf=imp:1", "gf=imp:2");
  EXPECT_TRUE(contains(error_of(replaced(two, "[levels]\nimp[0,0]=eps1", "[levels]\nimp[0,1]=eps1")), "[levels] imp[0,1]: an impurity level belongs to a diagonal element"));
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "imp[0,0]=eps1", "imp[0,0]=eps2")), "[levels] imp[0,0]: parameter \"eps2\" is not in params"));
  EXPECT_TRUE(contains(error_of(replaced(siam_qs, "imp[0,0]=eps1", "imp[1,1]=eps1")), "[levels] imp[1,1]: imp[1,1] is outside block"));
}

// ---- Against the shipped template ----

TEST(Template, SiamQs) {
  auto const dir = template_dir("SIAM/QS");
  auto const m   = load_manifest(dir);
  auto const r   = check_against_template(m, dir);
  EXPECT_FALSE(r.field);
  EXPECT_EQ(r.correlator_type.at("A_d-A_d"), 'd');
  EXPECT_EQ(r.correlator_type.at("self_d-A_d"), 'd');
  EXPECT_EQ(r.correlator_type.at("n_d-n_d"), 's');
  EXPECT_EQ(r.correlator_type.at("sigma_d-sigma_d"), 't');
}

TEST(Template, Mismatches) {
  EXPECT_EQ(check_error_of(siam_qs), "");
  EXPECT_TRUE(contains(check_error_of(replaced(siam_qs, "symtype=QS", "symtype=QSZ")), "[template] symtype: \"QSZ\", but data.in states \"QS\""));
  EXPECT_TRUE(contains(check_error_of(replaced(siam_qs, "channels=1", "channels=2")), "[template] channels: 2, but data.in states 1"));
  EXPECT_TRUE(contains(check_error_of(replaced(siam_qs, "hop0 A_d", "hop0 hop1 A_d")), "operator \"hop1\" is not in data.in"));
  EXPECT_TRUE(contains(check_error_of(replaced(siam_qs, "hop0 A_d", "A_d")), "operator \"hop0\" of data.in is missing"));
  EXPECT_TRUE(contains(check_error_of(replaced(siam_qs, "=n_d-n_d", "=n_d-A_d")), "\"n_d\" (s) and \"A_d\" (d) form no NRG spectrum"));
  EXPECT_EQ(check_error_of(replaced(siam_qs, "imp[0,0]=SigmaHartree", "imp[0,0]=n_d^2")), ""); // any operator of type s
  EXPECT_TRUE(contains(check_error_of(replaced(siam_qs, "imp[0,0]=SigmaHartree", "imp[0,0]=A_d")), "operator \"A_d\" is not of type s"));
  auto const three = replaced(replaced(siam_qs, "gf=imp:1", "gf=imp:3"), "coefficient_channels=imp[0,0]", "coefficient_channels=imp[0,0] imp[1,1] imp[2,2]");
  EXPECT_TRUE(contains(check_error_of(three), "[bath] coefficient_channels: 3 entries, but data.in has 1 channels"));
}

TEST(Template, SpectrumLists) {
  // s with g is a specs correlator; p only pairs with p
  auto const g = modified_template("SIAM/QS", "siam_qs_g", "sSigmaHartree", "gSigmaHartree");
  auto const sg = replaced(replaced(siam_qs, "chi_NN imp[0,0]=n_d-n_d", "chi_NN imp[0,0]=n_d-SigmaHartree"), "[hartree]\nimp[0,0]=SigmaHartree\n", "");
  auto const sg_sigma = replaced(sg, "F_l imp[0,0]=self_d-A_d\nF_r imp[0,0]=A_d-self_d\nI imp[0,0]=self_d-self_d\n", "");
  auto const r = check_against_template(parse_manifest(sg_sigma, "manifest"), g);
  EXPECT_EQ(r.correlator_type.at("n_d-SigmaHartree"), 's');
  auto const p = modified_template("SIAM/QS", "siam_qs_p", "sn_d^2", "pn_d^2");
  EXPECT_TRUE(contains(check_error_of(replaced(siam_qs, "=n_d-n_d", "=n_d-n_d^2"), p), "(s) and \"n_d^2\" (p) form no NRG spectrum"));
  EXPECT_EQ(check_error_of(replaced(siam_qs, "=n_d-n_d", "=n_d^2-n_d^2"), p), "");
}

TEST(Template, SpinSuffix) {
  // QSZ has a field: d spectra are named per spin, others never are
  auto const qsz = template_dir("SIAM/QSZ");
  auto const r   = check_against_template(parse_manifest(siam_qsz, "manifest"), qsz);
  EXPECT_TRUE(r.field);
  EXPECT_EQ(r.correlator_type.at("A_d-A_d-u"), 'd');
  EXPECT_TRUE(contains(check_error_of(replaced(siam_qsz, "=A_d-A_d-u", "=A_d-A_d"), qsz), "has a field: a d spectrum needs the spin suffix"));
  EXPECT_TRUE(contains(check_error_of(replaced(siam_qsz, "=n_d-n_d", "=n_d-n_d-u"), qsz), "a spin suffix -u/-d belongs only to d spectra"));
  // QS has none
  EXPECT_TRUE(contains(check_error_of(replaced(siam_qs, "=A_d-A_d", "=A_d-A_d-u")), "a spin suffix -u/-d belongs only to d spectra"));
}

TEST(Template, ShippedManifests) {
  for (auto const *name : {"SIAM/QS", "SIAM/QSZ", "SIAM/ISO", "2orb-UJ/QS", "Holstein/Nph=10/QS", "CLEAN/QS"}) {
    auto const dir = template_dir(name);
    EXPECT_NO_THROW(check_against_template(load_manifest(dir), dir)) << name;
  }
}
