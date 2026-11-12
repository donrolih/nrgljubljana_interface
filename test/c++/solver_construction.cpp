// Constructing the solver from a template's manifest (no NRG run)

#include <nrgljubljana_interface/solver_core.hpp>

#include <triqs/test_tools/gfs.hpp>

#include <iterator>
#include <set>
#include <sstream>
#include <string>

using namespace nrgljubljana_interface;

namespace {

  solver_core make_solver(std::string const &model, std::string const &symtype) {
    constr_params_t cp;
    cp.model   = model;
    cp.symtype = symtype;
    return solver_core{cp};
  }

  // The names of a list value, without their order (NRG matches list entries by name)
  std::set<std::string> names(std::string const &list) {
    std::istringstream ss(list);
    return {std::istream_iterator<std::string>{ss}, std::istream_iterator<std::string>{}};
  }

} // namespace

TEST(SolverConstruction, SiamQs) {
  auto const S = make_solver("SIAM", "QS");
  EXPECT_EQ(S.gf_struct, (gf_struct_t{{"imp", 1}}));
  EXPECT_EQ(S.Delta_struct, S.gf_struct);
  EXPECT_EQ(S.chi_struct, (gf_struct_t{{"imp", 1}}));
  EXPECT_EQ(names(S.constr_params.ops), names("n_d n_d^2 hop0 A_d self_d sigma_d SigmaHartree"));
  EXPECT_EQ(names(S.constr_params.specs), names("n_d-n_d"));
  EXPECT_EQ(names(S.constr_params.specd), names("A_d-A_d self_d-A_d A_d-self_d self_d-self_d"));
  EXPECT_EQ(names(S.constr_params.spect), names("sigma_d-sigma_d"));
  EXPECT_TRUE(names(S.constr_params.specq).empty());
  EXPECT_TRUE(names(S.constr_params.specot).empty());
  EXPECT_EQ(names(S.constr_params.params), names("U1 eps1"));
  EXPECT_FALSE(S.constr_params.polarized);
}

TEST(SolverConstruction, SiamQszIsPolarized) {
  auto const S = make_solver("SIAM", "QSZ");
  EXPECT_EQ(S.gf_struct, (gf_struct_t{{"up", 1}, {"dn", 1}}));
  EXPECT_EQ(S.chi_struct, (gf_struct_t{{"imp", 1}}));
  EXPECT_EQ(names(S.constr_params.specs), names("n_d-n_d SZd-SZd"));
  EXPECT_TRUE(S.constr_params.polarized);
}

TEST(SolverConstruction, CleanHasNoG) {
  auto const S = make_solver("CLEAN", "QS");
  EXPECT_TRUE(S.gf_struct.empty());
  EXPECT_EQ(S.Delta_struct, (gf_struct_t{{"imp", 1}}));
  EXPECT_TRUE(S.chi_struct.empty());
}

TEST(SolverConstruction, NoManifest) {
  try {
    make_solver("2orb-UJt", "QS");
    FAIL() << "a template without a manifest was accepted";
  } catch (triqs::runtime_error const &e) { EXPECT_NE(std::string(e.what()).find("has no manifest"), std::string::npos) << e.what(); }
}

MAKE_MAIN
