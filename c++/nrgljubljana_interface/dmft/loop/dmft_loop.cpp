#include "./dmft_loop.hpp"

#include "../self_energy.hpp"

#include <mpi/mpi.hpp>
#include <nda/mpi.hpp>

#include <algorithm>
#include <cmath>
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

  } // namespace

  dmft_loop::dmft_loop(impurity_solver_t &solver, lattice_t const &lattice, mixer_t &mixer, convergence_t &convergence, loop_params_t params,
                       std::optional<mu_control_t> const &mu_control)
     : solver_(solver),
       lattice_(lattice),
       mixer_(mixer),
       convergence_(convergence),
       params_(std::move(params)),
       mu_control_(mu_control.value_or(mu_control_t{})) {
    auto const &m = solver_.manifest();
    check_options(params_, m);
    check_manifest(m, params_.solve_params);
    if (not(std::isfinite(params_.gamma_floor) and params_.gamma_floor >= 0) or not(std::isfinite(params_.sigma_floor) and params_.sigma_floor >= 0))
      throw std::invalid_argument("dmft_loop: the floors must be finite and >= 0");
    onsite_    = lattice_.onsite_energy(m.gf);
    degeneracy_ = m.degeneracy;
    if (not params_.run_dir.empty() and solver_.communicator().rank() == 0) std::filesystem::create_directories(params_.run_dir);
  }

  g_w_t const &dmft_loop::Delta_used() const {
    if (not Delta_used_) throw std::logic_error("dmft_loop: not initialized");
    return *Delta_used_;
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

  void dmft_loop::initialize(g_w_t const &Sigma0, double mu0) {
    if (not std::isfinite(mu0)) throw std::invalid_argument("dmft_loop: μ0 must be finite, got " + number(mu0));
    auto const &m = solver_.manifest();
    g_w_t const reference{solver_.mesh(), m.gf};
    detail::check_same_structure(Sigma0, reference, "dmft_loop", "Sigma0", "G");
    auto const Sigma         = causal_projection(Sigma0, params_.sigma_floor).first;
    auto const [G0, Delta0]  = lattice_.self_consistency(Sigma, mu0);
    Delta_used_             = hybridization_from_gamma(spectral_part(Delta0), onsite_, params_.gamma_floor).first;
    mu_                      = mu0;
    iteration_               = 0;
    last_mu_step_            = 0;
    Sigma_                   = Sigma;
    G_loc_                   = G0;
    G_loc_previous_.reset();
    last_.reset();
  }

  step_result_t dmft_loop::step() {
    if (not Delta_used_) throw std::logic_error("dmft_loop: step before initialize");
    auto const comm = solver_.communicator();
    long const n    = iteration_ + 1;
    step_result_t r;
    r.iteration = n;
    r.mu_used   = mu_;

    // 1–3: the levels, the solve, a causal Σ
    auto const &result = solver_.solve(*Delta_used_, with_levels(mu_), params_.options);
    auto const found   = result.processed.find("Sigma");
    if (found == result.processed.end()) throw std::runtime_error("dmft_loop: the solve gave no Σ");
    r.expv                 = result.expv;
    auto [Sigma, sigma_rep] = causal_projection(found->second, params_.sigma_floor);
    r.sigma_projection     = sigma_rep;

    // 4–6: G_loc at μ_used, μ_next, Δ_raw
    auto const G_loc      = lattice_.local_gf(Sigma, mu_);
    auto const adjusted   = mu_control_.adjust({Sigma, mu_, params_.solve_params.T, degeneracy_}, lattice_);
    r.mu                  = adjusted.report;
    r.filling             = adjusted.filling;
    bool const mixer_mu   = mixer_.moves_mu(n);
    double const mu_raw   = mixer_mu ? mu_ : adjusted.mu_next;
    auto const Delta_raw  = lattice_.self_consistency(Sigma, mu_raw).second;

    // 7: convergence
    block_matrices_t h_imp;
    for (auto const &[block, size] : solver_.manifest().gf) {
      nda::matrix<std::complex<double>> h(size, size);
      h = 0;
      for (long i = 0; i < size; ++i) h(i, i) = -mu_;
      h_imp.push_back(std::move(h));
    }
    auto const G_imp     = impurity_gf(*Delta_used_, Sigma, h_imp);
    double const mu_step = mixer_mu ? last_mu_step_ : adjusted.mu_next - mu_;
    r.convergence        = convergence_.decide({n, G_loc, G_loc_previous_, G_imp, adjusted.filling, mu_control_.params().goal, mu_step});

    // 8: the next Δ, also after a terminal decision
    g_w_t Delta_next = *Delta_used_;
    double mu_next   = mu_;
    try {
      auto const mixed = mixer_.mix({n, spectral_part(*Delta_used_), spectral_part(Delta_raw), mu_, adjusted.filling});
      r.mix            = mixed.report;
      if (mixer_mu and not mixed.mu_next) throw std::runtime_error("dmft_loop: the mixer moves μ but gave none");
      mu_next             = mixer_mu ? *mixed.mu_next : adjusted.mu_next;
      auto [next, report] = hybridization_from_gamma(mixed.gamma_next, onsite_, params_.gamma_floor);
      Delta_next          = std::move(next);
      r.delta_projection  = report;
    } catch (std::exception const &) {
      if (r.convergence.decision == "continue") throw;
      Delta_next   = *Delta_used_;
      mu_next      = mu_;
      r.kept_delta = true;
    }

    // 9: rank 0's Δ_next, μ_next and decision on every rank, then commit
    if (comm.size() > 1) {
      for (long b = 0; b < long(Delta_next.size()); ++b) mpi::broadcast(Delta_next[b].data(), comm, 0);
      mpi::broadcast(mu_next, comm, 0);
      mpi::broadcast(r.convergence.decision, comm, 0);
    }
    r.mu_next       = mu_next;
    last_mu_step_   = mu_next - mu_;
    mu_             = mu_next;
    Delta_used_     = std::move(Delta_next);
    G_loc_previous_ = G_loc;
    G_loc_          = G_loc;
    Sigma_          = std::move(Sigma);
    iteration_      = n;
    last_           = r;

    if (comm.rank() == 0) write_statistics(r);
    if (callback_) callback_(*this, r);
    return r;
  }

  std::string dmft_loop::run(long max_steps) {
    if (max_steps < 0) throw std::invalid_argument("dmft_loop: max_steps must be >= 0");
    for (long taken = 0;;) {
      auto const r = step();
      ++taken;
      if (r.convergence.decision != "continue") return r.convergence.decision;
      if (max_steps > 0 and taken >= max_steps) return r.convergence.decision;
    }
  }

  void dmft_loop::write_statistics(step_result_t const &r) const {
    if (params_.run_dir.empty()) return;
    auto const path  = std::filesystem::path(params_.run_dir) / "stats.dat";
    bool const fresh = r.iteration == 1 or not std::filesystem::exists(path);
    std::ofstream out(path, fresh ? std::ios::trunc : std::ios::app);
    if (not out) throw std::runtime_error("dmft_loop: cannot write " + path.string());
    if (fresh)
      out << "# iteration mu_used mu_next filling distance window_max distance_imp_loc mix_residual sigma_projected delta_projected mu_status decision\n";
    out << std::setprecision(12) << r.iteration << ' ' << r.mu_used << ' ' << r.mu_next << ' ' << r.filling << ' '
        << (r.convergence.distance ? *r.convergence.distance : std::numeric_limits<double>::quiet_NaN()) << ' ' << r.convergence.window_max << ' '
        << r.convergence.distance_imp_loc << ' ' << r.mix.residual_norm << ' ' << r.sigma_projection.projected_points << ' '
        << r.delta_projection.projected_points << ' ' << r.mu.status << ' ' << r.convergence.decision << '\n';
  }

} // namespace nrgljubljana_interface::dmft
