// The DMFT loop (dmft/loop/dmft_loop) with a mock solver. No NRG run.
//
// The mock solver returns Σ = c Δ. On the Bethe lattice of half-bandwidth D, Δ = (D/2)² G_loc at self-consistency,
// and G_loc = g_D(ω + μ − c (D/2)² G_loc) is the Bethe Green's function of half-bandwidth D' = D √(1 + c): the fixed
// point is known in closed form. The loop rebuilds Δ by Kramers-Kronig on a coarse mesh, and the √ band edges of the
// Bethe lattice limit that: G_loc agrees with the closed form to 1e-2 in the integrated spectral distance; for the same
// reason the spectral weight on this mesh is 1 only to about 1e-2, and the μ control is given weight_tol = 1e-2.

#include <gtest/gtest.h>

#include <nrgljubljana_interface/dmft/lattice/bethe_lattice.hpp>
#include <nrgljubljana_interface/dmft/loop/dmft_loop.hpp>
#include <nrgljubljana_interface/dmft/measurements.hpp>
#include <nrgljubljana_interface/templates/manifest.hpp>

#include <mpi/mpi.hpp>
#include <nda/mpi.hpp>

#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace nrgljubljana_interface;
using namespace nrgljubljana_interface::dmft;
using triqs::mesh::refreq_log;
using dcomplex = std::complex<double>;

namespace {

  // MPI for the run on several ranks; nothing without an MPI launcher
  class mpi_environment : public ::testing::Environment {
   public:
    void SetUp() override {
      if (mpi::has_env and not mpi::is_initialized()) MPI_Init(nullptr, nullptr);
    }
    void TearDown() override {
      if (mpi::has_env and mpi::is_initialized() and not mpi::is_finalized()) MPI_Finalize();
    }
  };
  [[maybe_unused]] ::testing::Environment *const mpi_env = ::testing::AddGlobalTestEnvironment(new mpi_environment);

  std::string const siam_manifest = R"([template]
format=1
model=SIAM
symtype=QS
channels=1
ops=n_d A_d self_d SigmaHartree
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

[hartree]
imp[0,0]=SigmaHartree
)";

  std::string const levels_section = "\n[levels]\nimp[0,0]=eps1\n";

  refreq_log test_mesh() { return {1e-3, 5.0, 1.05}; }

  double const D = 1.0, c = 0.5;

  // Σ = c Δ; records the levels it was given
  class mock_solver : public impurity_solver_t {
   public:
    explicit mock_solver(std::string const &manifest_text, mpi::communicator comm = {})
       : manifest_(templates::parse_manifest(manifest_text, "mock")), comm_(comm) {}

    result_t const &solve(g_w_t const &Delta, solve_params_t const &sp, options_t const &options) override {
      result_         = result_t{};
      result_.options = options;
      result_.mesh    = mesh();
      g_w_t Sigma     = c * Delta;
      result_.processed.emplace("Sigma", std::move(Sigma));
      result_.expv["n_d"] = 0.5;
      levels.push_back(sp.model_parameters.at("eps1"));
      return result_;
    }
    [[nodiscard]] templates::manifest_t const &manifest() const override { return manifest_; }
    [[nodiscard]] refreq_log mesh() const override { return test_mesh(); }
    [[nodiscard]] mpi::communicator communicator() const override { return comm_; }

    std::vector<double> levels;

   private:
    templates::manifest_t manifest_;
    mpi::communicator comm_;
    result_t result_;
  };

  loop_params_t loop_params(std::string run_dir = {}) {
    loop_params_t p;
    p.solve_params.Nz                = 1;
    p.solve_params.T                 = 0.01;
    p.solve_params.model_parameters = {{"U1", 1.0}};
    p.run_dir                        = std::move(run_dir);
    return p;
  }

  convergence_params_t convergence_params() {
    convergence_params_t p;
    p.eps      = 1e-9;
    p.window   = 2;
    p.max_iter = 200;
    return p;
  }

  g_w_t constant(dcomplex value) {
    g_w_t S{test_mesh(), triqs::gfs::gf_struct_t{{"imp", 1}}};
    S[0].data() = value;
    return S;
  }

  // The closed form of the fixed point: the Bethe G of half-bandwidth D' at ω + μ
  g_w_t fixed_point(double mu) {
    double const Dp = D * std::sqrt(1 + c);
    g_w_t G         = constant(0);
    auto const mesh = test_mesh();
    auto const &pts = mesh.points();
    for (long k = 0; k < long(pts.size()); ++k) {
      dcomplex const z    = {pts[k] + mu, 1e-12};
      G[0].data()(k, 0, 0) = 2.0 / (z + std::sqrt(z - Dp) * std::sqrt(z + Dp));
    }
    return G;
  }

  double max_abs_diff(g_w_t const &a, g_w_t const &b) { return double(max_element(abs(a[0].data() - b[0].data()))); }

  mu_control_params_t goal_params(double goal) {
    mu_control_params_t p;
    p.goal       = goal;
    p.weight_tol = 1e-2;
    return p;
  }

  template <typename Exception, typename Function> void expect_error(Function f, std::string const &detail) {
    try {
      f();
    } catch (Exception const &e) {
      EXPECT_NE(std::string(e.what()).find(detail), std::string::npos) << "message: " << e.what();
      return;
    }
    ADD_FAILURE() << "no exception containing: " << detail;
  }

} // namespace

TEST(Loop, FixedMu) {
  mock_solver solver(siam_manifest + levels_section);
  bethe_lattice_t const lattice(D);
  linear_mixer mixer(0.5);
  convergence_t convergence(convergence_params());
  dmft_loop loop(solver, lattice, mixer, convergence, loop_params());
  double const mu0 = 0.1;
  loop.initialize(constant(dcomplex(0.2, -1e-3)), mu0);
  EXPECT_EQ(loop.run(), "converged");
  EXPECT_EQ(loop.mu(), mu0);
  // Integrated: pointwise the deviation peaks at the √ band edges (0.04 measured), where the closed form is steepest
  EXPECT_LT(spectral_distance(*loop.G_loc(), fixed_point(mu0)).max, 1e-2);
  ASSERT_FALSE(solver.levels.empty());
  for (double level : solver.levels) EXPECT_EQ(level, -mu0); // on-site energy 0
  EXPECT_EQ(long(solver.levels.size()), loop.iteration());
}

// With a filling goal the converged μ is the root of n_D'(μ) = goal, here found by bisection on the closed form
TEST(Loop, FillingGoal) {
  mock_solver solver(siam_manifest + levels_section);
  bethe_lattice_t const lattice(D);
  linear_mixer mixer(0.5);
  convergence_t convergence(convergence_params());
  double const goal = 0.8;
  dmft_loop loop(solver, lattice, mixer, convergence, loop_params(), mu_control_t{goal_params(goal)});
  loop.initialize(constant(dcomplex(0.2, -1e-3)), 0.0);
  EXPECT_EQ(loop.run(), "converged");
  EXPECT_NEAR(loop.last_step()->filling, goal, 1e-6);

  std::map<std::string, long> const deg{{"imp", 2}};
  auto filling = [&](double mu) { return occupancy(fixed_point(mu), 0.01, deg).total; };
  double lo = -1, hi = 1;
  for (int i = 0; i < 50; ++i) {
    double const mid = (lo + hi) / 2;
    (filling(mid) < goal ? lo : hi) = mid;
  }
  EXPECT_NEAR(loop.mu(), (lo + hi) / 2, 1e-2);
}

TEST(Loop, StepRunCallbackStatistics) {
  std::string const dir = "dmft_loop_run";
  std::filesystem::remove_all(dir);
  auto make = [&](mock_solver &solver, linear_mixer &mixer, convergence_t &convergence, bethe_lattice_t const &lattice, std::string run_dir) {
    dmft_loop loop(solver, lattice, mixer, convergence, loop_params(std::move(run_dir)), mu_control_t{goal_params(0.8)});
    loop.initialize(constant(dcomplex(0.2, -1e-3)), 0.0);
    return loop;
  };
  bethe_lattice_t const lattice(D);
  mock_solver s1(siam_manifest + levels_section), s2(siam_manifest + levels_section);
  linear_mixer m1(0.5), m2(0.5);
  convergence_t c1(convergence_params()), c2(convergence_params());

  auto by_run = make(s1, m1, c1, lattice, dir);
  long calls  = 0;
  by_run.set_callback([&](dmft_loop const &, step_result_t const &) { ++calls; });
  EXPECT_EQ(by_run.run(5), "continue");
  EXPECT_EQ(by_run.iteration(), 5);
  EXPECT_EQ(calls, 5);

  auto by_step = make(s2, m2, c2, lattice, {});
  for (int i = 0; i < 5; ++i) (void)by_step.step();
  EXPECT_EQ(by_step.mu(), by_run.mu());
  EXPECT_EQ(max_abs_diff(by_step.Delta(), by_run.Delta()), 0.0);

  std::vector<std::string> lines;
  {
    std::ifstream stats(std::filesystem::path(dir) / "stats.dat");
    for (std::string line; std::getline(stats, line);) lines.push_back(line);
  }
  ASSERT_EQ(lines.size(), 6u);
  EXPECT_EQ(lines[0].rfind("# iteration mu mu_next", 0), 0u);
  EXPECT_EQ(lines[5].rfind("5 ", 0), 0u);
  std::filesystem::remove_all(dir);
}

TEST(Loop, Errors) {
  bethe_lattice_t const lattice(D);
  linear_mixer mixer(0.5);
  convergence_t convergence(convergence_params());

  mock_solver no_levels(siam_manifest);
  expect_error<std::invalid_argument>([&] { (void)dmft_loop(no_levels, lattice, mixer, convergence, loop_params()); },
                                      "no [levels] entry for imp[0,0]");
  mock_solver solver(siam_manifest + levels_section);
  auto p                         = loop_params();
  p.solve_params.model_parameters["eps1"] = -0.3;
  expect_error<std::invalid_argument>([&] { (void)dmft_loop(solver, lattice, mixer, convergence, p); }, "the model parameter eps1 is an impurity level");
  p                 = loop_params();
  p.options.broaden = false;
  expect_error<std::invalid_argument>([&] { (void)dmft_loop(solver, lattice, mixer, convergence, p); }, "must have fdm and broaden");
  p                    = loop_params();
  p.solve_params.Nz    = 2;
  p.options.z          = {0.5};
  expect_error<std::invalid_argument>([&] { (void)dmft_loop(solver, lattice, mixer, convergence, p); }, "the full set of z");
  p                     = loop_params();
  p.options.correlators = std::vector<std::string>{"G"};
  expect_error<std::invalid_argument>([&] { (void)dmft_loop(solver, lattice, mixer, convergence, p); }, "must select F_l");

  dmft_loop loop(solver, lattice, mixer, convergence, loop_params());
  expect_error<std::logic_error>([&] { (void)loop.step(); }, "step before initialize");
  g_w_t wrong{test_mesh(), triqs::gfs::gf_struct_t{{"up", 1}}};
  expect_error<std::invalid_argument>([&] { loop.initialize(wrong, 0.0); }, "Sigma0 does not have the blocks of G");
}

// Every rank holds rank 0's state after each step, and the same as a run on one rank
TEST(Loop, AllRanksAgree) {
  mpi::communicator world;
  bethe_lattice_t const lattice(D);
  auto run = [&](mpi::communicator comm, long steps) {
    mock_solver solver(siam_manifest + levels_section, comm);
    linear_mixer mixer(0.5);
    convergence_t convergence(convergence_params());
    dmft_loop loop(solver, lattice, mixer, convergence, loop_params(), mu_control_t{goal_params(0.8)});
    loop.initialize(constant(dcomplex(0.2, -1e-3)), 0.0);
    std::vector<g_w_t> deltas;
    std::vector<double> mus;
    for (long i = 0; i < steps; ++i) {
      (void)loop.step();
      deltas.push_back(loop.Delta());
      mus.push_back(loop.mu());
    }
    return std::pair{deltas, mus};
  };
  auto const [together, mu_together] = run(world, 4);
  auto const [alone, mu_alone]       = run(mpi::communicator{MPI_COMM_SELF}, 4);
  for (std::size_t i = 0; i < together.size(); ++i) {
    auto rank0 = together[i][0].data();
    mpi::broadcast(rank0, world, 0);
    EXPECT_EQ(double(max_element(abs(together[i][0].data() - rank0))), 0.0) << "step " << i + 1 << " on rank " << world.rank();
    EXPECT_EQ(double(max_element(abs(together[i][0].data() - alone[i][0].data()))), 0.0) << "step " << i + 1;
    EXPECT_EQ(mu_together[i], mu_alone[i]);
  }
}
