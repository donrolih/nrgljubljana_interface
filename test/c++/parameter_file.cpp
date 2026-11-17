// Parser for parameter-file text (nrg/parameter_file.hpp)

#include <gtest/gtest.h>

#include <nrgljubljana_interface/nrg/parameter_file.hpp>

#include <stdexcept>
#include <string>
#include <vector>

using namespace nrgljubljana_interface::nrg;

namespace {

  // The error message of read_parameter_sections(text, source), or "" if it parses
  std::string error_of(std::string const &text, std::string const &source = {}) {
    try {
      read_parameter_sections(text, source);
    } catch (std::invalid_argument const &e) { return e.what(); }
    return {};
  }

  bool contains(std::string const &s, std::string const &part) { return s.find(part) != std::string::npos; }

} // namespace

// ---- Valid input ----

TEST(ParameterFile, SectionsAndKeys) {
  auto s = read_parameter_sections("[a]\nx=1\ny=two words\n[b]\nx=3\n");
  ASSERT_EQ(s.size(), 2u);
  EXPECT_EQ(s.at("a").at("x"), "1");
  EXPECT_EQ(s.at("a").at("y"), "two words");
  EXPECT_EQ(s.at("b").at("x"), "3"); // the same key in another section is a different key
}

TEST(ParameterFile, WhitespaceCommentsAndBlankLines) {
  auto s = read_parameter_sections("# leading comment\n\n  [ sec ]  \n\t key  =  value \t\n   # indented comment\n\nk2=v2\r\n");
  ASSERT_EQ(s.size(), 1u);
  auto const &sec = s.at("sec");
  EXPECT_EQ(sec.size(), 2u);
  EXPECT_EQ(sec.at("key"), "value");
  EXPECT_EQ(sec.at("k2"), "v2"); // Windows line ending
}

TEST(ParameterFile, SplitAtFirstEquals) {
  auto s = read_parameter_sections("[c]\nG imp[0,0]=d:A_d-A_d\nexpr=a=b\n");
  EXPECT_EQ(s.at("c").at("G imp[0,0]"), "d:A_d-A_d"); // keys may contain spaces and brackets
  EXPECT_EQ(s.at("c").at("expr"), "a=b");
}

TEST(ParameterFile, EmptyValuesAndSections) {
  auto s = read_parameter_sections("[a]\nempty=\nspaces=   \n[b]\n");
  EXPECT_EQ(s.at("a").at("empty"), "");
  EXPECT_EQ(s.at("a").at("spaces"), "");
  EXPECT_TRUE(s.at("b").empty());
  EXPECT_TRUE(read_parameter_sections("").empty());
  EXPECT_TRUE(read_parameter_sections("# only a comment\n\n").empty());
}

TEST(ParameterFile, HashAfterValueIsPartOfValue) {
  // Comments are whole lines only, as in NRG Ljubljana
  auto s = read_parameter_sections("[a]\nx=1 # not a comment\n");
  EXPECT_EQ(s.at("a").at("x"), "1 # not a comment");
}

TEST(ParameterFile, NoFinalNewline) {
  auto s = read_parameter_sections("[a]\nx=1");
  EXPECT_EQ(s.at("a").at("x"), "1");
}

TEST(ParameterFile, Manifest) {
  // A SIAM/QS template manifest
  auto s = read_parameter_sections(R"([template]
format=1
model=SIAM
symtype=QS
channels=1
ops=n_d n_d^2 hop0 A_d self_d sigma_d SigmaHartree
# model parameters required in [extra]
params=U1 eps1

[structure]
gf=imp:1
delta=imp:1
chi=imp:1
degeneracy=imp:2

[bath]
layout=uncoupled
coefficient_channels=imp[0,0]

[correlators]
# <correlator> <block element> = <type>:<NRG correlator>
G imp[0,0]=d:A_d-A_d
F_l imp[0,0]=d:self_d-A_d
F_r imp[0,0]=d:A_d-self_d
I imp[0,0]=d:self_d-self_d
chi_NN imp[0,0]=s:n_d-n_d
chi_SS imp[0,0]=t:sigma_d-sigma_d

[hartree]
imp[0,0]=SigmaHartree

[levels]
imp[0,0]=eps1
)",
                                   "templates/SIAM/QS/manifest");
  EXPECT_EQ(s.size(), 6u);
  EXPECT_EQ(s.at("template").at("symtype"), "QS");
  EXPECT_EQ(split_list(s.at("template").at("ops")).size(), 7u);
  EXPECT_EQ(s.at("correlators").size(), 6u);
  EXPECT_EQ(s.at("correlators").at("F_l imp[0,0]"), "d:self_d-A_d");
  EXPECT_EQ(s.at("hartree").at("imp[0,0]"), "SigmaHartree");
}

// ---- Lists ----

TEST(ParameterFile, SplitList) {
  EXPECT_EQ(split_list("U1 eps1"), (std::vector<std::string>{"U1", "eps1"}));
  EXPECT_EQ(split_list("  up:1 \t dn:1  "), (std::vector<std::string>{"up:1", "dn:1"}));
  EXPECT_EQ(split_list("one"), (std::vector<std::string>{"one"}));
  EXPECT_TRUE(split_list("").empty());
  EXPECT_TRUE(split_list(" \t ").empty());
}

// ---- Invalid input: each error names the line ----

TEST(ParameterFile, LineWithoutEquals) {
  auto e = error_of("[a]\nx=1\ntypo line\n");
  EXPECT_TRUE(contains(e, "line 3")) << e;
  EXPECT_TRUE(contains(e, "typo line")) << e;
}

TEST(ParameterFile, MalformedSectionHeader) {
  EXPECT_TRUE(contains(error_of("[a\nx=1\n"), "line 1: malformed section header"));
  EXPECT_TRUE(contains(error_of("[a] trailing\n"), "malformed section header"));
  EXPECT_TRUE(contains(error_of("[]\n"), "empty section name"));
  EXPECT_TRUE(contains(error_of("[  ]\n"), "empty section name"));
}

TEST(ParameterFile, DuplicateSection) {
  auto e = error_of("[a]\nx=1\n[b]\n[a]\ny=2\n");
  EXPECT_TRUE(contains(e, "line 4")) << e;
  EXPECT_TRUE(contains(e, "[a] given twice")) << e;
}

TEST(ParameterFile, KeyBeforeSection) {
  auto e = error_of("# comment\nx=1\n[a]\n");
  EXPECT_TRUE(contains(e, "line 2: key before the first section")) << e;
}

TEST(ParameterFile, EmptyKey) {
  auto e = error_of("[a]\n=value\n");
  EXPECT_TRUE(contains(e, "line 2: empty key")) << e;
  EXPECT_TRUE(contains(error_of("[a]\n  = value\n"), "empty key"));
}

TEST(ParameterFile, DuplicateKey) {
  auto e = error_of("[a]\nx=1\n x = 2\n");
  EXPECT_TRUE(contains(e, "line 3")) << e;
  EXPECT_TRUE(contains(e, "\"x\" given twice")) << e;
}

TEST(ParameterFile, SourceInMessage) {
  auto e = error_of("[a]\nbad\n", "templates/SIAM/QS/manifest");
  EXPECT_EQ(e.rfind("templates/SIAM/QS/manifest: line 2:", 0), 0u) << e;
  // Without a source the message starts with the line
  EXPECT_EQ(error_of("[a]\nbad\n").rfind("line 2:", 0), 0u);
}
