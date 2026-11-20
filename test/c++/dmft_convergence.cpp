// Convergence (dmft/loop/convergence). No NRG run.
//
// G_loc of iteration k is a fixed Green's function times (1 + c_k): the spectral distance between two iterations is
// |c_k − c_(k−1)| times that of the function itself, so a sequence of distances is set through the c_k.

#include <gtest/gtest.h>

#include <nrgljubljana_interface/dmft/loop/convergence.hpp>
#include <nrgljubljana_interface/dmft/measurements.hpp>

#include <complex>
#include <filesystem>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

using namespace nrgljubljana_interface;
using namespace nrgljubljana_interface::dmft;
using triqs::mesh::refreq_log;
using dcomplex = std::complex<double>;

namespace {

  refreq_log test_mesh() { return {1e-2, 5.0, 1.2}; }

  g_w_t base() {
    auto const mesh = test_mesh();
    g_w_t G{mesh, triqs::gfs::gf_struct_t{{"imp", 1}}};
    auto const &pts = mesh.points();
    for (long k = 0; k < long(pts.size()); ++k) G[0].data()(k, 0, 0) = 1.0 / (pts[k] - 0.2 + dcomplex(0, 0.5));
    return G;
  }

  g_w_t scaled(double c) {
    g_w_t G = base();
    G[0].data() *= (1 + c);
    return G;
  }

  // ∫ |A| of the base function: the distance per unit of c
  double unit_distance() { return spectral_distance(scaled(1.0), base()).max; }

  // Runs the iterations 1, 2, … with the given distances between consecutive G_loc (the first entry is ignored)
  struct run_t {
    convergence_t convergence;
    double c = 0;
    long iteration = 0;
    std::optional<g_w_t> previous;

    convergence_result_t step(double distance, double filling = 0.8, std::optional<double> goal = {}, double mu_step = 0,
                              double imp_offset = 0) {
      ++iteration;
      if (iteration > 1) c += distance / unit_distance();
      g_w_t const G = scaled(c);
      auto const r  = convergence.decide({iteration, G, previous, scaled(c + imp_offset), filling, goal, mu_step});
      previous      = G;
      return r;
    }
  };

  convergence_params_t params(long min_iter = 1, long max_iter = 50) {
    convergence_params_t p;
    p.eps      = 1e-6;
    p.window   = 3;
    p.min_iter = min_iter;
    p.max_iter = max_iter;
    return p;
  }

  template <typename Function> void expect_invalid_argument(Function f, std::string const &detail) {
    try {
      f();
    } catch (std::invalid_argument const &e) {
      EXPECT_NE(std::string(e.what()).find(detail), std::string::npos) << "message: " << e.what();
      return;
    }
    ADD_FAILURE() << "no std::invalid_argument containing: " << detail;
  }

} // namespace

TEST(Convergence, Window) {
  run_t run{convergence_t{params()}};
  auto const r1 = run.step(0);
  EXPECT_FALSE(r1.distance.has_value());
  EXPECT_EQ(r1.decision, "continue");
  EXPECT_EQ(r1.window_samples, 0);
  std::vector<std::string> decisions;
  for (double d : {1e-3, 5e-7, 5e-7, 5e-7, 2e-6, 5e-7, 5e-7, 5e-7}) decisions.push_back(run.step(d).decision);
  EXPECT_EQ(decisions, (std::vector<std::string>{"continue", "continue", "continue", "converged", "continue", "continue", "continue", "converged"}));
  EXPECT_NEAR(*run.step(5e-7).distance, 5e-7, 1e-12);
  EXPECT_EQ(run.convergence.history().size(), 9u);
}

TEST(Convergence, Iterations) {
  {
    run_t run{convergence_t{params(6)}};
    std::vector<std::string> decisions;
    for (double d : {0.0, 5e-7, 5e-7, 5e-7, 5e-7, 5e-7}) decisions.push_back(run.step(d).decision);
    EXPECT_EQ(decisions, (std::vector<std::string>{"continue", "continue", "continue", "continue", "continue", "converged"}));
  }
  {
    run_t run{convergence_t{params(1, 3)}};
    std::vector<std::string> decisions;
    for (double d : {0.0, 1e-3, 1e-3}) decisions.push_back(run.step(d).decision);
    EXPECT_EQ(decisions.back(), "stop");
  }
  {
    run_t run{convergence_t{params(1, 4)}};
    std::vector<std::string> decisions;
    for (double d : {0.0, 5e-7, 5e-7, 5e-7}) decisions.push_back(run.step(d).decision);
    EXPECT_EQ(decisions.back(), "converged"); // converging at max_iter is not a stop
  }
}

TEST(Convergence, FillingAndMu) {
  auto last = [](double filling, std::optional<double> goal, double mu_step) {
    run_t run{convergence_t{params()}};
    convergence_result_t r;
    for (double d : {0.0, 5e-7, 5e-7, 5e-7}) r = run.step(d, filling, goal, mu_step);
    return r;
  };
  EXPECT_EQ(last(0.8, 0.8, 0.0).decision, "converged");
  auto const off = last(0.80001, 0.8, 0.0);
  EXPECT_EQ(off.decision, "continue");
  EXPECT_FALSE(off.filling_ok);
  auto const moving = last(0.8, 0.8, 1e-6);
  EXPECT_EQ(moving.decision, "continue");
  EXPECT_FALSE(moving.mu_ok);
  EXPECT_EQ(last(0.80001, std::nullopt, 1e-6).decision, "converged"); // fixed μ: not checked
}

TEST(Convergence, ImpurityLatticeGate) {
  auto last = [](std::optional<double> gate) {
    auto p        = params();
    p.eps_imp_loc = gate;
    run_t run{convergence_t{p}};
    convergence_result_t r;
    double const offset = 1e-4 / unit_distance();
    for (double d : {0.0, 5e-7, 5e-7, 5e-7}) r = run.step(d, 0.8, std::nullopt, 0.0, offset);
    return r;
  };
  auto const ungated = last(std::nullopt);
  EXPECT_EQ(ungated.decision, "converged");
  EXPECT_NEAR(ungated.distance_imp_loc, 1e-4, 1e-12);
  auto const gated = last(1e-5);
  EXPECT_EQ(gated.decision, "continue");
  EXPECT_FALSE(gated.imp_loc_ok);
}

// Interrupted after three iterations, written and read into a new object: the same decisions as without
TEST(Convergence, Checkpoint) {
  std::vector<double> const distances{0.0, 1e-3, 5e-7, 5e-7, 5e-7, 5e-7};
  run_t full{convergence_t{params()}};
  std::vector<std::string> expected;
  for (double d : distances) expected.push_back(full.step(d).decision);

  std::string const path = "dmft_convergence.h5";
  run_t first{convergence_t{params()}};
  std::vector<std::string> decisions;
  for (std::size_t i = 0; i < 3; ++i) decisions.push_back(first.step(distances[i]).decision);
  {
    h5::file f(path, 'w');
    first.convergence.h5_write(h5::group(f), "convergence");
  }
  run_t second{convergence_t{params()}, first.c, first.iteration, first.previous};
  {
    h5::file f(path, 'r');
    second.convergence.h5_read(h5::group(f), "convergence");
    convergence_params_t p;
    h5_read(h5::group(f).open_group("convergence"), "params", p);
    EXPECT_EQ(p.window, 3);
  }
  for (std::size_t i = 3; i < distances.size(); ++i) decisions.push_back(second.step(distances[i]).decision);
  EXPECT_EQ(decisions, expected);
  std::filesystem::remove(path);
}

TEST(Convergence, Errors) {
  g_w_t const G = base();
  convergence_t c{params()};
  expect_invalid_argument([&] { c.decide({0, G, std::nullopt, G, 0.8, std::nullopt, 0}); }, "the iteration must be at least 1");
  expect_invalid_argument([&] { c.decide({2, G, std::nullopt, G, 0.8, std::nullopt, 0}); }, "iteration 2 needs the G_loc of the previous iteration");
  (void)c.decide({1, G, std::nullopt, G, 0.8, std::nullopt, 0});
  expect_invalid_argument([&] { c.decide({3, G, G, G, 0.8, std::nullopt, 0}); }, "iteration 3 does not follow iteration 1");
  g_w_t other{test_mesh(), triqs::gfs::gf_struct_t{{"up", 1}}};
  expect_invalid_argument([&] { c.decide({2, G, G, other, 0.8, std::nullopt, 0}); }, "G_imp does not have the blocks of G_loc");

  auto bad = [](auto change, std::string const &detail) {
    auto p = params();
    change(p);
    expect_invalid_argument([&] { convergence_t{p}; }, detail);
  };
  bad([](auto &p) { p.eps = 0; }, "eps must be positive");
  bad([](auto &p) { p.window = 0; }, "window must be at least 1");
  bad([](auto &p) { p.max_iter = 0; }, "need 1 ≤ min_iter ≤ max_iter");
  bad([](auto &p) { p.mu_tol = std::numeric_limits<double>::quiet_NaN(); }, "mu_tol must be positive");
  bad([](auto &p) { p.eps_imp_loc = -1.0; }, "eps_imp_loc must be positive");
}
