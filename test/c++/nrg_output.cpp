// Output files of NRG Ljubljana through nrg/output.hpp (NRG Ljubljana's readers, the interface's types); the fixtures
// are taken from a SIAM/QS run

#include <gtest/gtest.h>

#include <nrgljubljana_interface/nrg/output.hpp>

#include <stdexcept>
#include <string>
#include <vector>

using namespace nrgljubljana_interface::nrg;

namespace {

  std::string const custom = R"(#               1                2                3                4                5
#               T     SigmaHartree             hop0              n_d            n_d^2
      2.164042561     0.2494574998    -0.1450298996     0.9978299992       1.46701563
      1.082021281     0.2489896891    -0.2819629262     0.9959587566      1.435962978
)";

  std::string const tdfdm = R"(#               T            E_fdm            C_fdm            F_fdm            S_fdm
            0.001     -2.807840478      2.071807846     -2.814620886      6.780408095
)";

  std::string const annotated = R"(0 (0 1) [1]
0.13862699 (-1 2) [2]
0.14154254 (1 2) [2]

0 (0 2) [2]
0.55414531 (-1 3) (1 3) [6]

)";

  template <typename F> std::string error_of(F f) {
    try {
      f();
    } catch (std::runtime_error const &e) { return e.what(); }
    return {};
  }

  bool contains(std::string const &s, std::string const &part) { return s.find(part) != std::string::npos; }

} // namespace

TEST(NrgOutput, StepTable) {
  auto const t = parse_step_table(custom);
  EXPECT_EQ(t.columns, (std::vector<std::string>{"T", "SigmaHartree", "hop0", "n_d", "n_d^2"})); // not the numbers
  ASSERT_EQ(t.values.extent(0), 2);
  ASSERT_EQ(t.values.extent(1), 5);
  EXPECT_EQ(t.values(0, 0), 2.164042561);
  EXPECT_EQ(t.values(1, 3), 0.9959587566);
  EXPECT_EQ(t.values(1, 2), -0.2819629262);
}

TEST(NrgOutput, FinalValues) {
  auto const v = parse_final_values(tdfdm);
  EXPECT_EQ(v.size(), 5u);
  EXPECT_EQ(v.at("T"), 0.001);
  EXPECT_EQ(v.at("F_fdm"), -2.814620886);
  EXPECT_TRUE(contains(error_of([] { parse_final_values(custom, "customfdm"); }), "customfdm: expected exactly one row, found 2"));
}

TEST(NrgOutput, StepTableErrors) {
  EXPECT_TRUE(contains(error_of([] { parse_step_table("1 2\n", "td"); }), "td: line 1: data before a header line"));
  EXPECT_TRUE(contains(error_of([] { parse_step_table("# T E\n1 2 3\n", "td"); }), "td: line 2: expected 2 values, found 3"));
  EXPECT_TRUE(contains(error_of([] { parse_step_table("# T E\n1 x\n", "td"); }), "not a number: \"x\""));
  EXPECT_TRUE(contains(error_of([] { parse_step_table("", "td"); }), "no header line"));
  EXPECT_TRUE(contains(error_of([] { parse_step_table("# T E\n1 2\n# T E\n", "td"); }), "header line after the data"));
}

TEST(NrgOutput, EnergyFlow) {
  auto const flow = parse_energy_flow(annotated);
  ASSERT_EQ(flow.size(), 2u);
  ASSERT_EQ(flow[0].size(), 3u);
  EXPECT_EQ(flow[0][0], (energy_level_t{0.0, {{0, 1}}, 1}));
  EXPECT_EQ(flow[0][1], (energy_level_t{0.13862699, {{-1, 2}}, 2}));
  ASSERT_EQ(flow[1].size(), 2u);
  EXPECT_EQ(flow[1][1], (energy_level_t{0.55414531, {{-1, 3}, {1, 3}}, 6}));
}

TEST(NrgOutput, EnergyFlowEmptyStepsAndEnd) {
  // A step without levels is a lone empty line; a missing final empty line still ends the last step
  auto const flow = parse_energy_flow("0 (0 1) [1]\n\n\n0 (0 2) [2]");
  ASSERT_EQ(flow.size(), 3u);
  EXPECT_EQ(flow[0].size(), 1u);
  EXPECT_TRUE(flow[1].empty());
  EXPECT_EQ(flow[2].size(), 1u);
}

TEST(NrgOutput, EnergyFlowErrors) {
  for (auto const *bad : {"x (0 1) [1]\n", "0 [1]\n", "0 (0 1\n", "0 (0 a) [1]\n", "0 (0 1)\n", "0 (0 1) [0]\n", "0 (0 1) [1] x\n"})
    EXPECT_TRUE(contains(error_of([&] { parse_energy_flow(bad, "annotated.dat"); }), "annotated.dat: line 1: ")) << bad;
}

TEST(NrgOutput, SpectrumFile) {
  EXPECT_EQ(spectrum_file("A_d", "A_d", "", 'd'), "spec_FDM_dens_A_d-A_d.bin");
  EXPECT_EQ(spectrum_file("A_d", "A_d", "u", 'd'), "spec_FDM_dens_A_d-A_d-u.bin");
  EXPECT_EQ(spectrum_file("A_d", "A_d", "d", 'd'), "spec_FDM_dens_A_d-A_d-d.bin");
  EXPECT_EQ(spectrum_file("n_d", "n_d", "", 's'), "corr_FDM_dens_n_d-n_d.bin");
  EXPECT_EQ(spectrum_file("sigma_d", "sigma_d", "", 't'), "spin_FDM_dens_sigma_d-sigma_d.bin");
  EXPECT_EQ(spectrum_file("x", "x", "", 'o'), "orbspin_FDM_dens_x-x.bin");
  EXPECT_EQ(spectrum_file("x", "x", "", 'q'), "specq_FDM_dens_x-x.bin");
  // operator names may contain '-': the split is the manifest's, not guessed from the name
  EXPECT_EQ(spectrum_file("n_d", "SigmaHartree-u", "", 's'), "corr_FDM_dens_n_d-SigmaHartree-u.bin");
  EXPECT_THROW(spectrum_file("x", "x", "", 'p'), std::invalid_argument);
  EXPECT_THROW(spectrum_file("A_d", "A_d", "x", 'd'), std::invalid_argument);
  EXPECT_THROW(spectrum_file("", "A_d", "", 'd'), std::invalid_argument);
}
