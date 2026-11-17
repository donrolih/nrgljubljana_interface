// The NRG parameter file of a solve and how the options set it (solver/param_file.hpp)

#include <gtest/gtest.h>

#include <nrgljubljana_interface/nrg/parameter_file.hpp>
#include <nrgljubljana_interface/solver/param_file.hpp>
#include <nrgljubljana_interface/templates/manifest.hpp>

#include <cstdlib>
#include <iterator>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace nrgljubljana_interface;

namespace {

  struct template_t {
    templates::manifest_t manifest;
    std::map<std::string, char> correlator_type;
  };

  template_t load(std::string const &name) {
    auto const *root = std::getenv("NRGIF_TEMPLATE_DIR");
    if (root == nullptr) throw std::runtime_error("NRGIF_TEMPLATE_DIR is not set");
    auto const dir = std::string(root) + "/" + name;
    auto m         = templates::load_manifest(dir);
    auto check     = templates::check_against_template(m, dir);
    return {m, check.correlator_type};
  }

  // The names of a list value, without their order (NRG matches list entries by name)
  std::set<std::string> names(std::string const &list) {
    std::istringstream ss(list);
    return {std::istream_iterator<std::string>{ss}, std::istream_iterator<std::string>{}};
  }

  struct rendered_t {
    constr_params_t cp;
    nrg_params_t np;
    nrg::parameter_sections_t file; // the rendered parameter file, parsed
  };

  rendered_t render(template_t const &t, options_t const &r, nrg_params_t np = {}) {
    constr_params_t cp;
    cp.model   = t.manifest.model;
    cp.symtype = t.manifest.symtype;
    solve_params_t sp;
    sp.model_parameters = {{"U1", 0.5}, {"eps1", -0.25}};
    apply_options(r, t.manifest, t.correlator_type, cp, np);
    return {cp, np, nrg::read_parameter_sections(render_param_file(cp, sp, np), "param")};
  }

} // namespace

TEST(ParamRendering, DefaultOptions) {
  auto const t = load("SIAM/QS");
  auto const r = render(t, {});
  auto const &p = r.file.at("param");
  EXPECT_EQ(names(p.at("specd")), names("A_d-A_d self_d-A_d A_d-self_d self_d-self_d"));
  EXPECT_EQ(names(p.at("specs")), names("n_d-n_d"));
  EXPECT_EQ(names(p.at("spect")), names("sigma_d-sigma_d"));
  EXPECT_EQ(p.at("specq"), "");
  EXPECT_EQ(p.at("specot"), "");
  EXPECT_EQ(p.at("fdm"), "true");
  EXPECT_EQ(p.at("fdmexpv"), "true");
  EXPECT_EQ(p.at("dmnrg"), "false");
  EXPECT_EQ(p.at("cfs"), "false");
  EXPECT_EQ(p.at("savebins"), "true");
  EXPECT_EQ(p.at("dumpannotated"), "0");
  EXPECT_EQ(p.at("model"), "SIAM");
  EXPECT_EQ(p.at("symtype"), "QS");
  EXPECT_EQ(r.file.at("extra").at("U1"), "0.5");
  EXPECT_EQ(r.file.at("extra").at("eps1"), "-0.25");
}

TEST(ParamRendering, ForwardOnly) {
  auto const t = load("SIAM/QS");
  nrg_params_t np;
  np.dmnrg = true;
  np.cfs   = true;
  options_t req;
  req.fdm      = false;
  req.energies = 20;
  auto const p = render(t, req, np).file.at("param");
  for (auto const *key : {"fdm", "fdmexpv", "dmnrg", "cfs", "savebins"}) EXPECT_EQ(p.at(key), "false") << key;
  for (auto const *key : {"specs", "specd", "spect", "specq", "specot"}) EXPECT_EQ(p.at(key), "") << key;
  EXPECT_EQ(p.at("dumpannotated"), "20");
}

TEST(ParamRendering, SubsetOfCorrelators) {
  auto const t = load("SIAM/QS");
  options_t req;
  req.correlators = std::vector<std::string>{"G", "chi_NN"};
  auto const p    = render(t, req).file.at("param");
  EXPECT_EQ(p.at("specd"), "A_d-A_d");
  EXPECT_EQ(p.at("specs"), "n_d-n_d");
  EXPECT_EQ(p.at("spect"), "");
  EXPECT_EQ(p.at("savebins"), "true");
}

TEST(ParamRendering, FdmWithoutCorrelators) {
  auto const t = load("SIAM/QS");
  nrg_params_t np;
  np.fdmexpv = false; // the options decide
  options_t req;
  req.correlators = std::vector<std::string>{};
  auto const p    = render(t, req, np).file.at("param");
  EXPECT_EQ(p.at("fdm"), "true");
  EXPECT_EQ(p.at("fdmexpv"), "true");
  EXPECT_EQ(p.at("savebins"), "false");
  EXPECT_EQ(p.at("specd"), "");
}

TEST(ParamRendering, EnergiesZeroKeepsLowLevelValue) {
  auto const t = load("SIAM/QS");
  nrg_params_t np;
  np.dumpannotated = 5;
  EXPECT_EQ(render(t, {}, np).file.at("param").at("dumpannotated"), "5");
}

TEST(ParamRendering, PolarizedSpectra) {
  auto const t = load("SIAM/QSZ");
  options_t req;
  req.correlators = std::vector<std::string>{"G", "chi_SS"};
  auto const p    = render(t, req).file.at("param");
  EXPECT_EQ(names(p.at("specd")), names("A_d-A_d-u A_d-A_d-d"));
  EXPECT_EQ(p.at("specs"), "SZd-SZd");
}

TEST(ParamRendering, TwistValues) {
  EXPECT_EQ(twist_values({}, 4), (std::vector<double>{0.25, 0.5, 0.75, 1.0}));
  EXPECT_EQ(twist_values({.z = {1.0}}, 4), (std::vector<double>{1.0}));
  EXPECT_TRUE(complete_z_average(twist_values({}, 4), 4));
  EXPECT_TRUE(complete_z_average({1.0, 0.5}, 2));
  EXPECT_FALSE(complete_z_average({1.0}, 2));
  EXPECT_FALSE(complete_z_average({0.25, 1.0}, 2));
  EXPECT_TRUE(complete_z_average({1.0}, 1));
}

TEST(ParamRendering, Method) {
  solve_params_t sp;
  EXPECT_FALSE(check_method(sp, {}));
  sp.method = "fdm";
  EXPECT_TRUE(check_method(sp, {}));
  EXPECT_THROW(check_method(sp, {.fdm = false}), std::invalid_argument);
  sp.method = "dmnrg";
  EXPECT_THROW(check_method(sp, {}), std::invalid_argument);
}
