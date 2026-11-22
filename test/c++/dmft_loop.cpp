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
#include <nrgljubljana_interface/solver/result.hpp>
#include <nrgljubljana_interface/templates/manifest.hpp>

#include <mpi/mpi.hpp>
#include <nda/mpi.hpp>

#include <cmath>
#include <cstdio>
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
      twist_result_t twist;
      twist.spectra["A_d-A_d"] = nrg::spectrum_bins_t{{-0.5, 0.5}, {0.5, 0.5}, {}, false};
      result_.twists.push_back(std::move(twist));
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

// Broyden reaches the fixed point of linear mixing in fewer iterations
TEST(Loop, Broyden) {
  bethe_lattice_t const lattice(D);
  double const mu0 = 0.1;
  auto iterations = [&](mixer_t &mixer) {
    mock_solver solver(siam_manifest + levels_section);
    convergence_t convergence(convergence_params());
    dmft_loop loop(solver, lattice, mixer, convergence, loop_params());
    loop.initialize(constant(dcomplex(0.2, -1e-3)), mu0);
    EXPECT_EQ(loop.run(), "converged");
    EXPECT_LT(spectral_distance(*loop.G_loc(), fixed_point(mu0)).max, 1e-2);
    return loop.iteration();
  };
  linear_mixer linear(0.5);
  broyden_mixer broyden;
  long const n_linear = iterations(linear), n_broyden = iterations(broyden);
  EXPECT_LT(n_broyden, n_linear) << "Broyden " << n_broyden << ", linear " << n_linear;
}

// With mu, the μ control moves μ during the warmup and only measures the filling afterwards
TEST(Loop, BroydenMovesMu) {
  mock_solver solver(siam_manifest + levels_section);
  bethe_lattice_t const lattice(D);
  broyden_params_t p;
  p.mu     = true;
  p.warmup = 2;
  broyden_mixer mixer(p);
  convergence_t convergence(convergence_params());
  double const goal = 0.8;
  dmft_loop loop(solver, lattice, mixer, convergence, loop_params(), mu_control_t{goal_params(goal)});
  std::vector<std::string> statuses;
  loop.set_callback([&](step_result_t const &r) { statuses.push_back(r.mu_report.status); });
  loop.initialize(constant(dcomplex(0.2, -1e-3)), 0.0);
  EXPECT_EQ(loop.run(), "converged");
  EXPECT_NEAR(loop.last_step()->filling, goal, 1e-6);
  ASSERT_GT(statuses.size(), 2u);
  EXPECT_NE(statuses[0], "mixer");
  EXPECT_NE(statuses[1], "mixer");
  for (std::size_t i = 2; i < statuses.size(); ++i) EXPECT_EQ(statuses[i], "mixer");

  broyden_mixer other(p);
  convergence_t c2(convergence_params());
  expect_error<std::invalid_argument>([&] { (void)dmft_loop(solver, lattice, other, c2, loop_params()); },
                                      "the mixer moves μ, which needs a μ control with a goal");
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
  by_run.set_callback([&](step_result_t const &) { ++calls; });
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

namespace {

  // On several ranks rank 0 prepares the directory, the others wait
  void fresh_directory(std::string const &dir, mpi::communicator const &world) {
    if (world.rank() == 0) std::filesystem::remove_all(dir);
    world.barrier();
  }

  struct run_record_t {
    std::vector<double> mus;
    std::vector<std::string> statuses;
  };

  // A loop with a filling goal in run_dir; the controls are owned here so a new object stands for a new process
  struct loop_set_t {
    mock_solver solver;
    bethe_lattice_t lattice{D};
    linear_mixer mixer;
    convergence_t convergence{convergence_params()};
    dmft_loop loop;
    run_record_t record;

    loop_set_t(std::string const &run_dir, double alpha = 0.5, double T = 0.01, bool keep_history = true)
       : solver(siam_manifest + levels_section), mixer(alpha), loop(make(run_dir, T, keep_history)) {
      loop.set_callback([this](step_result_t const &r) {
        record.mus.push_back(r.mu_next);
        record.statuses.push_back(r.convergence.status);
      });
    }

    dmft_loop make(std::string const &run_dir, double T, bool keep_history) {
      auto p            = loop_params(run_dir);
      p.solve_params.T  = T;
      p.keep_history    = keep_history;
      return dmft_loop(solver, lattice, mixer, convergence, p, mu_control_t{goal_params(0.8)});
    }
  };

  result_t read_result(std::string const &run_dir) {
    h5::file f((std::filesystem::path(run_dir) / "checkpoint.h5").string(), 'r');
    result_t r;
    h5_read(h5::group(f), "result", r);
    return r;
  }

} // namespace

// Interrupted after 3 iterations and resumed with new objects: the same iterations as an uninterrupted run
TEST(Loop, Resume) {
  mpi::communicator world;
  std::string const full_dir = "dmft_loop_full", cut_dir = "dmft_loop_cut";
  fresh_directory(full_dir, world);
  fresh_directory(cut_dir, world);

  loop_set_t full(full_dir);
  full.loop.initialize(constant(dcomplex(0.2, -1e-3)), 0.0);
  EXPECT_EQ(full.loop.run(), "converged");
  long const n = full.loop.iteration();
  ASSERT_GT(n, 4);

  run_record_t cut_record;
  {
    loop_set_t first(cut_dir);
    first.loop.initialize(constant(dcomplex(0.2, -1e-3)), 0.0);
    EXPECT_EQ(first.loop.run(3), "continue");
    cut_record = first.record;
    EXPECT_TRUE(read_result(cut_dir).twists.at(0).spectra.empty()); // no bins before the end
  }
  if (world.rank() == 0) std::ofstream(std::filesystem::path(cut_dir) / "checkpoint.h5.tmp") << "left over from an interrupted write";
  world.barrier();

  loop_set_t second(cut_dir);
  second.loop.resume(cut_dir);
  EXPECT_EQ(second.loop.iteration(), 3);
  EXPECT_EQ(second.loop.run(), "converged");
  EXPECT_EQ(second.loop.iteration(), n);
  auto mus = cut_record.mus;
  mus.insert(mus.end(), second.record.mus.begin(), second.record.mus.end());
  EXPECT_EQ(mus, full.record.mus);
  EXPECT_EQ(max_abs_diff(second.loop.Delta(), full.loop.Delta()), 0.0);

  // Bins only in the terminal checkpoint
  EXPECT_FALSE(read_result(full_dir).twists.at(0).spectra.empty());

  // One history entry per iteration; the last G_loc is the loop's
  {
    h5::file f((std::filesystem::path(full_dir) / "history.h5").string(), 'r');
    h5::group const root(f);
    EXPECT_TRUE(root.has_subgroup("iteration_0001"));
    char last[32];
    std::snprintf(last, sizeof last, "iteration_%04ld", n);
    ASSERT_TRUE(root.has_subgroup(last));
    g_w_t G;
    h5::read(root.open_group(last), "G_loc", G);
    EXPECT_EQ(max_abs_diff(G, *full.loop.G_loc()), 0.0);
  }

  // A new run from the converged one starts at its Δ_mix and μ_next, with fresh controls
  std::string const next_dir = "dmft_loop_next";
  fresh_directory(next_dir, world);
  loop_set_t next(next_dir);
  next.loop.initialize_from(full_dir);
  EXPECT_EQ(next.loop.iteration(), 0);
  EXPECT_EQ(next.loop.mu(), full.loop.mu());
  EXPECT_EQ(max_abs_diff(next.loop.Delta(), full.loop.Delta()), 0.0);
  EXPECT_TRUE(next.convergence.history().empty());

  world.barrier();
  if (world.rank() == 0)
    for (auto const &dir : {full_dir, cut_dir, next_dir}) std::filesystem::remove_all(dir);
}

TEST(Loop, ResumeChecksParameters) {
  mpi::communicator world;
  std::string const dir = "dmft_loop_check";
  fresh_directory(dir, world);
  {
    loop_set_t first(dir);
    first.loop.initialize(constant(dcomplex(0.2, -1e-3)), 0.0);
    (void)first.loop.run(2);
  }
  {
    loop_set_t other(dir, 0.4, 0.02);
    expect_error<std::invalid_argument>([&] { other.loop.resume(dir); }, "the solve parameters; the mixer (stored \"linear alpha=0.5\"");
    other.loop.resume(dir, true);
    EXPECT_EQ(other.loop.iteration(), 2);
  }
  {
    loop_set_t empty("dmft_loop_none");
    expect_error<std::invalid_argument>([&] { empty.loop.resume("dmft_loop_none"); }, "no checkpoint.h5 in dmft_loop_none");
  }
  {
    std::string const quiet = "dmft_loop_quiet";
    fresh_directory(quiet, world);
    loop_set_t no_history(quiet, 0.5, 0.01, false);
    no_history.loop.initialize(constant(dcomplex(0.2, -1e-3)), 0.0);
    (void)no_history.loop.run(2);
    EXPECT_TRUE(std::filesystem::exists(std::filesystem::path(quiet) / "checkpoint.h5"));
    EXPECT_FALSE(std::filesystem::exists(std::filesystem::path(quiet) / "history.h5"));
    world.barrier();
    if (world.rank() == 0) std::filesystem::remove_all(quiet);
  }
  world.barrier();
  if (world.rank() == 0) {
    std::filesystem::remove_all(dir);
    std::filesystem::remove_all("dmft_loop_none");
  }
}

