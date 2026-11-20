#include "./dmft_loop.hpp"

#include "../self_energy.hpp"

#include <mpi/mpi.hpp>
#include <nda/mpi.hpp>

#include <fcntl.h>
#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

namespace nrgljubljana_interface::dmft {

  namespace {

    using detail::number;

    long block_index(templates::block_structure_t const &s, std::string const &name) {
      for (long b = 0; b < long(s.size()); ++b)
        if (s[b].first == name) return b;
      throw std::invalid_argument("dmft_loop: block " + name + " is not in the structure of G");
    }

    void check_options(loop_params_t const &p, templates::manifest_t const &m) {
      std::string const pre = "dmft_loop: ";
      validate_options(p.options, m);
      if (not p.options.fdm or not p.options.broaden) throw std::invalid_argument(pre + "the options must have fdm and broaden");
      auto const selected = selected_correlators(p.options, m);
      for (auto const *name : templates::gf_correlators)
        if (std::find(selected.begin(), selected.end(), name) == selected.end())
          throw std::invalid_argument(pre + "the options must select " + name + " (G, F_l, F_r and I give Σ)");
      if (not p.options.z.empty()) {
        auto z = p.options.z;
        std::sort(z.begin(), z.end());
        bool full = long(z.size()) == long(p.solve_params.Nz);
        for (std::size_t i = 0; full and i < z.size(); ++i) full = std::abs(z[i] - double(i + 1) / p.solve_params.Nz) <= 1e-12;
        if (not full) throw std::invalid_argument(pre + "the options must solve the full set of z, i/Nz for i = 1..Nz");
      }
    }

    void check_manifest(templates::manifest_t const &m, solve_params_t const &sp) {
      std::string const pre = "dmft_loop: ";
      if (not m.has_self_energy()) throw std::invalid_argument(pre + "the template does not give G, F_l, F_r and I");
      if (m.delta != m.gf) throw std::invalid_argument(pre + "the template's Δ and G have different block structures");
      for (auto const &[block, n] : m.gf)
        for (long i = 0; i < n; ++i)
          if (not m.levels.contains(templates::block_element_t{block, i, i}))
            throw std::invalid_argument(pre + "the manifest has no [levels] entry for " + block + "[" + std::to_string(i) + "," + std::to_string(i) + "]");
      for (auto const &[element, name] : m.levels)
        if (sp.model_parameters.contains(name))
          throw std::invalid_argument(pre + "the model parameter " + name + " is an impurity level, set by the loop from the on-site energy and μ");
    }

    void write_reports(h5::group g, step_result_t const &r) {
      h5_write(g, "Sigma_projection", r.Sigma_projection);
      h5_write(g, "Delta_projection", r.Delta_projection);
      h5_write(g, "mu", r.mu_report);
      h5_write(g, "mix", r.mix_report);
      h5_write(g, "convergence", r.convergence);
      h5::write(g, "expv", r.expv);
    }

    void read_reports(h5::group g, step_result_t &r) {
      h5_read(g, "Sigma_projection", r.Sigma_projection);
      h5_read(g, "Delta_projection", r.Delta_projection);
      h5_read(g, "mu", r.mu_report);
      h5_read(g, "mix", r.mix_report);
      h5_read(g, "convergence", r.convergence);
      r.expv.clear();
      h5::read(g, "expv", r.expv);
    }

    // Run f on rank 0 only; an exception there is raised on every rank, so no rank is left waiting in a collective
    template <typename F> void on_rank_0(mpi::communicator const &comm, F &&f) {
      std::string error;
      if (comm.rank() == 0) {
        try {
          f();
        } catch (std::exception const &e) { error = e.what(); }
      }
      if (comm.size() > 1) mpi::broadcast(error, comm, 0);
      if (not error.empty()) throw std::runtime_error(error);
    }

    // Data on disk before the rename that publishes the file
    void sync(std::filesystem::path const &path) {
      int const fd = ::open(path.c_str(), O_RDONLY);
      if (fd < 0) throw std::runtime_error("dmft_loop: cannot open " + path.string() + " to flush it");
      ::fsync(fd);
      ::close(fd);
    }

    std::filesystem::path checkpoint_path(std::string const &run_dir) { return std::filesystem::path(run_dir) / "checkpoint.h5"; }

    std::string history_name(long iteration) {
      char name[32];
      std::snprintf(name, sizeof name, "iteration_%04ld", iteration);
      return name;
    }

    // The checkpoint of a run, or std::invalid_argument
    h5::file open_checkpoint(std::string const &run_dir) {
      auto const path = checkpoint_path(run_dir);
      if (not std::filesystem::exists(path)) throw std::invalid_argument("dmft_loop: no checkpoint.h5 in " + run_dir);
      return h5::file(path.string(), 'r');
    }

  } // namespace

  dmft_loop::dmft_loop(impurity_solver_t &solver, lattice_t const &lattice, mixer_t &mixer, convergence_t &convergence, loop_params_t params,
                       std::optional<mu_control_t> const &mu_control)
     : solver_(solver),
       lattice_(lattice),
       mixer_(mixer),
       convergence_(convergence),
       params_(std::move(params)),
       mu_control_(mu_control.value_or(mu_control_t{})) {
    setup();
  }

  dmft_loop::dmft_loop(solver_core &solver, lattice_t const &lattice, mixer_t &mixer, convergence_t &convergence, loop_params_t params,
                       std::optional<mu_control_t> const &mu_control)
     : owned_solver_(std::make_unique<solver_core_adapter>(solver)),
       solver_(*owned_solver_),
       lattice_(lattice),
       mixer_(mixer),
       convergence_(convergence),
       params_(std::move(params)),
       mu_control_(mu_control.value_or(mu_control_t{})) {
    setup();
  }

  void dmft_loop::setup() {
    auto const &m = solver_.manifest();
    check_options(params_, m);
    check_manifest(m, params_.solve_params);
    if (not(std::isfinite(params_.Gamma_floor) and params_.Gamma_floor >= 0) or not(std::isfinite(params_.Sigma_floor) and params_.Sigma_floor >= 0))
      throw std::invalid_argument("dmft_loop: the floors must be finite and >= 0");
    onsite_    = lattice_.onsite_energy(m.gf);
    degeneracy_ = m.degeneracy;
    if (not params_.run_dir.empty()) on_rank_0(solver_.communicator(), [&] { std::filesystem::create_directories(params_.run_dir); });
  }

  g_w_t const &dmft_loop::Delta() const {
    if (not Delta_) throw std::logic_error("dmft_loop: not initialized");
    return *Delta_;
  }

  solve_params_t dmft_loop::with_levels(double mu) const {
    auto const &m = solver_.manifest();
    auto sp       = params_.solve_params;
    std::map<std::string, double> set;
    for (auto const &[e, name] : m.levels) {
      long const b       = block_index(m.gf, e.block);
      double const value = onsite_[b](e.i, e.i).real() - mu;
      if (auto const it = set.find(name); it != set.end()) {
        if (std::abs(it->second - value) > 1e-12 * std::max(1.0, std::abs(value)))
          throw std::invalid_argument("dmft_loop: the levels mapped to " + name + " differ: " + number(it->second) + " and " + number(value));
        continue;
      }
      set[name]                  = value;
      sp.model_parameters[name] = value;
    }
    return sp;
  }

  void dmft_loop::initialize(std::complex<double> Sigma0, double mu0) {
    g_w_t S{solver_.mesh(), solver_.manifest().gf};
    for (long b = 0; b < long(S.size()); ++b) {
      auto &d = S[b].data();
      d       = 0;
      for (long k = 0; k < d.extent(0); ++k)
        for (long i = 0; i < d.extent(1); ++i) d(k, i, i) = Sigma0;
    }
    initialize(S, mu0);
  }

  void dmft_loop::initialize(g_w_t const &Sigma0, double mu0) {
    if (not std::isfinite(mu0)) throw std::invalid_argument("dmft_loop: μ0 must be finite, got " + number(mu0));
    auto const &m = solver_.manifest();
    g_w_t const reference{solver_.mesh(), m.gf};
    detail::check_same_structure(Sigma0, reference, "dmft_loop", "Sigma0", "G");
    auto const Sigma         = causal_projection(Sigma0, params_.Sigma_floor).first;
    auto const [G0, Delta0]  = lattice_.self_consistency(Sigma, mu0);
    Delta_             = hybridization_from_Gamma(spectral_part(Delta0), onsite_, params_.Gamma_floor).first;
    mu_                      = mu0;
    iteration_               = 0;
    last_mu_step_            = 0;
    Sigma_                   = Sigma;
    G_loc_                   = G0;
    G_loc_previous_.reset();
    last_.reset();
    mixer_.reset();
    convergence_.reset();
    start_files();
  }

  void dmft_loop::start_files() const {
    if (params_.run_dir.empty()) return;
    on_rank_0(solver_.communicator(), [&] { std::filesystem::remove(std::filesystem::path(params_.run_dir) / "history.h5"); });
  }

  void dmft_loop::resume(std::string const &run_dir, bool allow_changes) {
    auto file = open_checkpoint(run_dir);
    h5::group const root(file);
    auto const params = root.open_group("params");
    std::vector<std::string> differences;
    auto stored_string = [&](std::string const &key) {
      std::string v;
      h5::read(params, key, v);
      return v;
    };
    solve_params_t sp;
    h5::read(params, "solve_params", sp);
    if (not(sp == params_.solve_params)) differences.emplace_back("the solve parameters");
    options_t options;
    h5::read(params, "options", options);
    if (not(options == params_.options)) differences.emplace_back("the options");
    double Gamma_floor = 0, Sigma_floor = 0;
    h5::read(params, "Gamma_floor", Gamma_floor);
    h5::read(params, "Sigma_floor", Sigma_floor);
    if (Gamma_floor != params_.Gamma_floor or Sigma_floor != params_.Sigma_floor) differences.emplace_back("the floors");
    mu_control_params_t mp;
    h5_read(params, "mu_control", mp);
    if (not(mp == mu_control_.params())) differences.emplace_back("the μ-control parameters");
    convergence_params_t cp;
    h5_read(params, "convergence", cp);
    if (not(cp == convergence_.params())) differences.emplace_back("the convergence parameters");
    if (auto const m = stored_string("mixer"); m != mixer_.description())
      differences.emplace_back("the mixer (stored \"" + m + "\", now \"" + mixer_.description() + "\")");
    if (auto const l = stored_string("lattice"); l != lattice_.description())
      differences.emplace_back("the lattice (stored \"" + l + "\", now \"" + lattice_.description() + "\")");
    if (stored_string("manifest") != solver_.manifest().text) differences.emplace_back("the template's manifest");
    triqs::mesh::refreq_log mesh;
    h5::read(params, "mesh", mesh);
    if (not(mesh == solver_.mesh())) differences.emplace_back("the mesh");
    if (not differences.empty() and not allow_changes) {
      std::string list;
      for (auto const &d : differences) list += (list.empty() ? "" : "; ") + d;
      throw std::invalid_argument("dmft_loop: the run in " + run_dir + " differs from this loop in " + list + " (allow_changes to continue anyway)");
    }

    auto const loop  = root.open_group("loop");
    auto const state = root.open_group("state");
    long iteration   = 0;
    double mu_next = 0, last_step = 0;
    h5::read(loop, "iteration", iteration);
    h5::read(loop, "mu_next", mu_next);
    h5::read(loop, "last_mu_step", last_step);
    g_w_t Delta_mix, Sigma, G_loc;
    h5::read(state, "Delta_mix", Delta_mix);
    h5::read(state, "Sigma", Sigma);
    h5::read(state, "G_loc", G_loc);
    auto const controls = root.open_group("controls");
    mixer_.h5_read(controls, "mixer");
    convergence_.h5_read(controls, "convergence");

    iteration_      = iteration;
    mu_             = mu_next;
    last_mu_step_   = last_step;
    Delta_          = std::move(Delta_mix);
    Sigma_          = std::move(Sigma);
    G_loc_          = G_loc;
    G_loc_previous_ = std::move(G_loc);
    last_.reset();
  }

  void dmft_loop::initialize_from(std::string const &run_dir) {
    auto file = open_checkpoint(run_dir);
    h5::group const root(file);
    auto const params = root.open_group("params");
    std::string manifest;
    h5::read(params, "manifest", manifest);
    triqs::mesh::refreq_log mesh;
    h5::read(params, "mesh", mesh);
    if (manifest != solver_.manifest().text or not(mesh == solver_.mesh()))
      throw std::invalid_argument("dmft_loop: the run in " + run_dir + " has another template or mesh");
    double mu_next = 0;
    h5::read(root.open_group("loop"), "mu_next", mu_next);
    g_w_t Delta_mix;
    h5::read(root.open_group("state"), "Delta_mix", Delta_mix);

    Delta_        = std::move(Delta_mix);
    mu_           = mu_next;
    iteration_    = 0;
    last_mu_step_ = 0;
    Sigma_.reset();
    G_loc_.reset();
    G_loc_previous_.reset();
    last_.reset();
    mixer_.reset();
    convergence_.reset();
    start_files();
  }

  step_result_t dmft_loop::step() {
    if (not Delta_) throw std::logic_error("dmft_loop: step before initialize");
    auto const comm = solver_.communicator();
    long const n    = iteration_ + 1;
    step_result_t r;
    r.iteration = n;
    r.mu        = mu_;

    // 1–3: the levels, the solve, a causal Σ
    auto const &result = solver_.solve(*Delta_, with_levels(mu_), params_.options);
    auto const found   = result.processed.find("Sigma");
    if (found == result.processed.end()) throw std::runtime_error("dmft_loop: the solve gave no Σ");
    r.expv                 = result.expv;
    auto [Sigma, Sigma_report] = causal_projection(found->second, params_.Sigma_floor);
    r.Sigma_projection     = Sigma_report;

    // 4–6: G_loc at μ, μ_next, Δ_out
    auto const G_loc      = lattice_.local_gf(Sigma, mu_);
    auto const adjusted   = mu_control_.adjust({Sigma, mu_, params_.solve_params.T, degeneracy_}, lattice_);
    r.mu_report           = adjusted.report;
    r.filling             = adjusted.filling;
    bool const mixer_mu   = mixer_.moves_mu(n);
    double const mu_out   = mixer_mu ? mu_ : adjusted.mu_next;
    auto const Delta_out  = lattice_.self_consistency(Sigma, mu_out).second;

    // 7: convergence
    block_matrices_t h_imp;
    for (auto const &[block, size] : solver_.manifest().gf) {
      nda::matrix<std::complex<double>> h(size, size);
      h = 0;
      for (long i = 0; i < size; ++i) h(i, i) = -mu_;
      h_imp.push_back(std::move(h));
    }
    auto const G_imp     = impurity_gf(*Delta_, Sigma, h_imp);
    double const mu_step = mixer_mu ? last_mu_step_ : adjusted.mu_next - mu_;
    r.convergence        = convergence_.decide({n, G_loc, G_loc_previous_, G_imp, adjusted.filling, mu_control_.params().goal, mu_step});

    // 8: the next Δ, also after a terminal status
    g_w_t Delta_mix = *Delta_;
    double mu_next   = mu_;
    try {
      auto const mixed = mixer_.mix({n, spectral_part(*Delta_), spectral_part(Delta_out), mu_, adjusted.filling});
      r.mix_report     = mixed.report;
      if (mixer_mu and not mixed.mu_next) throw std::runtime_error("dmft_loop: the mixer moves μ but gave none");
      mu_next             = mixer_mu ? *mixed.mu_next : adjusted.mu_next;
      auto [next, report] = hybridization_from_Gamma(mixed.Gamma_mix, onsite_, params_.Gamma_floor);
      Delta_mix          = std::move(next);
      r.Delta_projection  = report;
    } catch (std::exception const &) {
      if (r.convergence.status == "continue") throw;
      Delta_mix   = *Delta_;
      mu_next      = mu_;
      r.kept_Delta = true;
    }

    // 9: rank 0's Δ_mix, μ_next and status on every rank, then commit
    if (comm.size() > 1) {
      for (long b = 0; b < long(Delta_mix.size()); ++b) mpi::broadcast(Delta_mix[b].data(), comm, 0);
      mpi::broadcast(mu_next, comm, 0);
      mpi::broadcast(r.convergence.status, comm, 0);
    }
    r.mu_next = mu_next;
    if (not params_.run_dir.empty())
      on_rank_0(comm, [&] {
        write_checkpoint(r, result, Delta_mix, Sigma, G_loc, G_imp);
        if (params_.keep_history) append_history(r, Delta_out, Sigma, G_loc);
      });
    last_mu_step_   = mu_next - mu_;
    mu_             = mu_next;
    Delta_     = std::move(Delta_mix);
    G_loc_previous_ = G_loc;
    G_loc_          = G_loc;
    Sigma_          = std::move(Sigma);
    iteration_      = n;
    last_           = r;

    on_rank_0(comm, [&] { write_statistics(r); });
    if (callback_) callback_(r);
    return r;
  }

  std::string dmft_loop::run(long max_steps) {
    if (max_steps < 0) throw std::invalid_argument("dmft_loop: max_steps must be >= 0");
    for (long taken = 0;;) {
      auto const r = step();
      ++taken;
      if (r.convergence.status != "continue") return r.convergence.status;
      if (max_steps > 0 and taken >= max_steps) return r.convergence.status;
    }
  }

  void dmft_loop::write_statistics(step_result_t const &r) const {
    if (params_.run_dir.empty()) return;
    auto const path  = std::filesystem::path(params_.run_dir) / "stats.dat";
    bool const fresh = r.iteration == 1 or not std::filesystem::exists(path);
    std::ofstream out(path, fresh ? std::ios::trunc : std::ios::app);
    if (not out) throw std::runtime_error("dmft_loop: cannot write " + path.string());
    if (fresh)
      out << "# iteration mu mu_next filling distance window_max distance_imp_loc mix_residual Sigma_projected Delta_projected mu_status status\n";
    out << std::setprecision(12) << r.iteration << ' ' << r.mu << ' ' << r.mu_next << ' ' << r.filling << ' '
        << (r.convergence.distance ? *r.convergence.distance : std::numeric_limits<double>::quiet_NaN()) << ' ' << r.convergence.window_max << ' '
        << r.convergence.distance_imp_loc << ' ' << r.mix_report.residual_norm << ' ' << r.Sigma_projection.projected_points << ' '
        << r.Delta_projection.projected_points << ' ' << r.mu_report.status << ' ' << r.convergence.status << '\n';
  }

  void dmft_loop::write_checkpoint(step_result_t const &r, result_t const &result, g_w_t const &Delta_mix, g_w_t const &Sigma, g_w_t const &G_loc,
                                   g_w_t const &G_imp) const {
    auto const final_path = checkpoint_path(params_.run_dir);
    auto const tmp        = std::filesystem::path(final_path.string() + ".tmp");
    {
      h5::file f(tmp.string(), 'w');
      h5::group const root(f);
      auto loop = root.create_group("loop");
      h5::write(loop, "iteration", r.iteration);
      h5::write(loop, "mu", r.mu);
      h5::write(loop, "mu_next", r.mu_next);
      h5::write(loop, "last_mu_step", r.mu_next - r.mu);
      h5::write(loop, "filling", r.filling);
      h5::write(loop, "status", r.convergence.status);
      h5::write(loop, "kept_Delta", long(r.kept_Delta));

      auto params = root.create_group("params");
      h5::write(params, "solve_params", params_.solve_params);
      h5::write(params, "options", params_.options);
      h5::write(params, "Gamma_floor", params_.Gamma_floor);
      h5::write(params, "Sigma_floor", params_.Sigma_floor);
      h5_write(params, "mu_control", mu_control_.params());
      h5_write(params, "convergence", convergence_.params());
      h5::write(params, "mixer", mixer_.description());
      h5::write(params, "lattice", lattice_.description());
      h5::write(params, "manifest", solver_.manifest().text);
      h5::write(params, "mesh", solver_.mesh());

      auto state = root.create_group("state");
      h5::write(state, "Delta", *Delta_);
      h5::write(state, "Delta_mix", Delta_mix);
      h5::write(state, "Sigma", Sigma);
      h5::write(state, "G_loc", G_loc);
      h5::write(state, "G_imp", G_imp);

      auto controls = root.create_group("controls");
      mixer_.h5_write(controls, "mixer");
      convergence_.h5_write(controls, "convergence");

      write_reports(root.create_group("reports"), r);
      nrgljubljana_interface::h5_write(root, "result", result, r.convergence.status != "continue");
    }
    sync(tmp);
    std::filesystem::rename(tmp, final_path);
  }

  void dmft_loop::append_history(step_result_t const &r, g_w_t const &Delta_out, g_w_t const &Sigma, g_w_t const &G_loc) const {
    h5::file f((std::filesystem::path(params_.run_dir) / "history.h5").string(), 'a');
    auto g = h5::group(f).create_group(history_name(r.iteration));
    h5::write(g, "mu", r.mu);
    h5::write(g, "mu_next", r.mu_next);
    h5::write(g, "filling", r.filling);
    h5::write(g, "status", r.convergence.status);
    h5::write(g, "Delta", *Delta_);
    h5::write(g, "Delta_out", Delta_out);
    h5::write(g, "Sigma", Sigma);
    h5::write(g, "G_loc", G_loc);
    write_reports(g.create_group("reports"), r);
  }

  void h5_write(h5::group g, std::string const &name, loop_params_t const &p) {
    auto grp = g.create_group(name);
    h5::write_hdf5_format(grp, p);
    h5_write(grp, "solve_params", p.solve_params);
    h5_write(grp, "options", p.options);
    h5::write(grp, "Gamma_floor", p.Gamma_floor);
    h5::write(grp, "Sigma_floor", p.Sigma_floor);
    h5::write(grp, "run_dir", p.run_dir);
    h5::write(grp, "keep_history", long(p.keep_history));
  }

  void h5_read(h5::group g, std::string const &name, loop_params_t &p) {
    auto grp          = g.open_group(name);
    long keep_history = 0;
    h5_read(grp, "solve_params", p.solve_params);
    h5_read(grp, "options", p.options);
    h5::read(grp, "Gamma_floor", p.Gamma_floor);
    h5::read(grp, "Sigma_floor", p.Sigma_floor);
    h5::read(grp, "run_dir", p.run_dir);
    h5::read(grp, "keep_history", keep_history);
    p.keep_history = keep_history != 0;
  }

  void h5_write(h5::group g, std::string const &name, step_result_t const &r) {
    auto grp = g.create_group(name);
    h5::write_hdf5_format(grp, r);
    h5::write(grp, "iteration", r.iteration);
    h5::write(grp, "mu", r.mu);
    h5::write(grp, "mu_next", r.mu_next);
    h5::write(grp, "filling", r.filling);
    h5::write(grp, "kept_Delta", long(r.kept_Delta));
    write_reports(grp.create_group("reports"), r);
  }

  void h5_read(h5::group g, std::string const &name, step_result_t &r) {
    auto grp        = g.open_group(name);
    long kept_Delta = 0;
    h5::read(grp, "iteration", r.iteration);
    h5::read(grp, "mu", r.mu);
    h5::read(grp, "mu_next", r.mu_next);
    h5::read(grp, "filling", r.filling);
    h5::read(grp, "kept_Delta", kept_Delta);
    r.kept_Delta = kept_Delta != 0;
    read_reports(grp.open_group("reports"), r);
  }

} // namespace nrgljubljana_interface::dmft
