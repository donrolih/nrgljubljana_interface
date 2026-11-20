#include "./mu_control.hpp"

#include "../measurements.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <memory>
#include <stdexcept>

namespace nrgljubljana_interface::dmft {

  namespace {

    using detail::number;

    // Restores GSL's error handler on every exit
    struct gsl_handler_off {
      gsl_error_handler_t *previous = gsl_set_error_handler_off();
      ~gsl_handler_off() { gsl_set_error_handler(previous); }
      gsl_handler_off()                                   = default;
      gsl_handler_off(gsl_handler_off const &)            = delete;
      gsl_handler_off &operator=(gsl_handler_off const &) = delete;
    };

    struct evaluation_t {
      double filling, weight_error;
    };

  } // namespace

  void h5_write(h5::group g, std::string const &name, mu_control_params_t const &p) {
    auto grp = g.create_group(name);
    if (p.goal) h5::write(grp, "goal", *p.goal);
    h5::write(grp, "maxdx", p.maxdx);
    h5::write(grp, "under", p.under);
    h5::write(grp, "root_tol", p.root_tol);
    h5::write(grp, "mu_tol", p.mu_tol);
    h5::write(grp, "max_evaluations", p.max_evaluations);
    h5::write(grp, "weight_tol", p.weight_tol);
    h5::write(grp, "epsabs", p.epsabs);
    h5::write(grp, "epsrel", p.epsrel);
  }

  void h5_read(h5::group g, std::string const &name, mu_control_params_t &p) {
    auto grp = g.open_group(name);
    p.goal   = {};
    if (grp.has_key("goal")) {
      double goal = 0;
      h5::read(grp, "goal", goal);
      p.goal = goal;
    }
    h5::read(grp, "maxdx", p.maxdx);
    h5::read(grp, "under", p.under);
    h5::read(grp, "root_tol", p.root_tol);
    h5::read(grp, "mu_tol", p.mu_tol);
    h5::read(grp, "max_evaluations", p.max_evaluations);
    h5::read(grp, "weight_tol", p.weight_tol);
    h5::read(grp, "epsabs", p.epsabs);
    h5::read(grp, "epsrel", p.epsrel);
  }

  mu_control_t::mu_control_t(mu_control_params_t const &params) : params_(params) {
    auto const &p   = params_;
    auto positive   = [](double x) { return std::isfinite(x) and x > 0; };
    std::string pre = "mu_control_t: ";
    if (p.goal and not(std::isfinite(*p.goal) and *p.goal >= 0)) throw std::invalid_argument(pre + "the goal must be finite and >= 0, got " + number(*p.goal));
    if (not positive(p.maxdx)) throw std::invalid_argument(pre + "maxdx must be positive, got " + number(p.maxdx));
    if (not(p.under > 0 and p.under <= 1)) throw std::invalid_argument(pre + "under must satisfy 0 < under ≤ 1, got " + number(p.under));
    if (not positive(p.root_tol)) throw std::invalid_argument(pre + "root_tol must be positive, got " + number(p.root_tol));
    if (not positive(p.mu_tol)) throw std::invalid_argument(pre + "mu_tol must be positive, got " + number(p.mu_tol));
    if (p.max_evaluations < 4) throw std::invalid_argument(pre + "max_evaluations must be at least 4, got " + std::to_string(p.max_evaluations));
    if (not positive(p.weight_tol)) throw std::invalid_argument(pre + "weight_tol must be positive, got " + number(p.weight_tol));
    if (not(std::isfinite(p.epsabs) and p.epsabs >= 0 and std::isfinite(p.epsrel) and p.epsrel >= 0))
      throw std::invalid_argument(pre + "epsabs and epsrel must be finite and >= 0");
  }

  mu_result_t mu_control_t::adjust(mu_input_t const &input, lattice_t const &lattice) const {
    auto const &p = params_;
    if (not std::isfinite(input.mu_used)) throw std::invalid_argument("mu_control_t: μ_used must be finite, got " + number(input.mu_used));

    // The number of states, for the plausibility of the filling
    double states = 0;
    for (long b = 0; b < long(input.Sigma.size()); ++b) {
      auto const it = input.degeneracy.find(input.Sigma.block_names()[b]);
      if (it != input.degeneracy.end()) states += double(it->second) * double(input.Sigma[b].data().extent(1));
    }

    std::map<double, evaluation_t> cache;
    auto evaluate = [&](double mu, long reserve) -> evaluation_t const & {
      if (not std::isfinite(mu)) throw std::runtime_error("mu_control_t: a trial μ is not finite");
      if (auto const found = cache.find(mu); found != cache.end()) return found->second;
      if (long(cache.size()) >= p.max_evaluations - reserve)
        throw std::runtime_error("mu_control_t: max_evaluations (" + std::to_string(p.max_evaluations) + ") exceeded");
      auto const G      = lattice.local_gf(input.Sigma, mu);
      double weight_err = 0;
      auto const w      = spectral_weight(G);
      for (auto const &m : w)
        for (long i = 0; i < m.extent(0); ++i) weight_err = std::max(weight_err, std::abs(m(i, i).real() - 1));
      if (p.goal and not(weight_err <= p.weight_tol))
        throw std::runtime_error("mu_control_t: the lattice spectral weight at μ = " + number(mu) + " is off by " + number(weight_err)
                                 + ", more than weight_tol; the mesh may not cover the spectrum");
      double const n   = occupancy(G, input.T, input.degeneracy, p.epsabs, p.epsrel).total;
      double const tol = p.weight_tol * states;
      if (p.goal and not(n >= -tol and n <= states + tol)) throw std::runtime_error("mu_control_t: unphysical filling " + number(n) + " at μ = " + number(mu));
      return cache.emplace(mu, evaluation_t{n, weight_err}).first->second;
    };

    double const mu       = input.mu_used;
    auto const initial    = evaluate(mu, 0);
    mu_result_t result;
    result.filling                  = initial.filling;
    result.report.filling_used      = initial.filling;
    result.report.weight_error_used = initial.weight_error;

    auto finish = [&](double mu_next, std::string status) {
      auto const last                 = evaluate(mu_next, 0);
      result.mu_next                  = mu_next;
      result.report.status            = std::move(status);
      result.report.filling_next      = last.filling;
      result.report.weight_error_next = last.weight_error;
      result.report.evaluations      = long(cache.size());
      return result;
    };

    if (not p.goal) return finish(mu, "fixed");
    double const goal      = *p.goal;
    double const error_old = initial.filling - goal;
    if (std::abs(error_old) <= p.root_tol) return finish(mu, "within_tolerance");

    double const direction = error_old < 0 ? 1.0 : -1.0;
    double const edge_mu   = mu + direction * p.maxdx / p.under;
    if (not std::isfinite(edge_mu) or edge_mu == mu) throw std::runtime_error("mu_control_t: the bracket cannot be represented at μ = " + number(mu));
    auto const edge         = evaluate(edge_mu, 0);
    double const error_edge = edge.filling - goal;
    if (direction * (edge.filling - initial.filling) < -p.root_tol)
      throw std::runtime_error("mu_control_t: the filling did not increase with μ between " + number(mu) + " and " + number(edge_mu));

    double mu_next = 0;
    std::string status;
    if (error_edge == 0 or std::signbit(error_old) != std::signbit(error_edge)) {
      gsl_handler_off const handler;
      struct callback_t {
        decltype(evaluate) &eval;
        double goal;
        std::exception_ptr error;
      } callback{evaluate, goal, {}};
      gsl_function function;
      function.params   = &callback;
      function.function = [](double x, void *context) noexcept -> double {
        auto &c = *static_cast<callback_t *>(context);
        if (not c.error) {
          try {
            return c.eval(x, 1).filling - c.goal; // one evaluation kept for μ_next
          } catch (...) { c.error = std::current_exception(); }
        }
        return std::numeric_limits<double>::quiet_NaN();
      };
      auto checked = [&](int code) {
        if (callback.error) std::rethrow_exception(callback.error);
        if (code != GSL_SUCCESS) throw std::runtime_error(std::string("mu_control_t: Brent's method failed: ") + gsl_strerror(code));
      };
      std::unique_ptr<gsl_root_fsolver, decltype(&gsl_root_fsolver_free)> solver(gsl_root_fsolver_alloc(gsl_root_fsolver_brent),
                                                                                  &gsl_root_fsolver_free);
      if (not solver) throw std::bad_alloc();
      checked(gsl_root_fsolver_set(solver.get(), &function, std::min(mu, edge_mu), std::max(mu, edge_mu)));

      double root    = edge_mu;
      bool converged = error_edge == 0;
      long const iterations = p.max_evaluations - long(cache.size()) - 1;
      for (long it = 0; not converged and it < iterations; ++it) {
        checked(gsl_root_fsolver_iterate(solver.get()));
        root                  = gsl_root_fsolver_root(solver.get());
        double const residual = evaluate(root, 1).filling - goal;
        converged             = residual == 0
           or gsl_root_test_interval(gsl_root_fsolver_x_lower(solver.get()), gsl_root_fsolver_x_upper(solver.get()), p.mu_tol,
                                     4 * std::numeric_limits<double>::epsilon())
              == GSL_SUCCESS;
      }
      if (not converged or not std::isfinite(root))
        throw std::runtime_error("mu_control_t: Brent's method did not converge within max_evaluations (" + std::to_string(p.max_evaluations) + ")");
      if (std::abs(evaluate(root, 1).filling - goal) > p.root_tol)
        throw std::runtime_error("mu_control_t: the filling at the root differs from the goal by more than root_tol");
      result.report.mu_root = root;
      mu_next               = mu + std::clamp(p.under * (root - mu), -p.maxdx, p.maxdx);
      status                = "bracketed";
    } else {
      mu_next = mu + direction * p.maxdx;
      status  = "limited";
    }
    // Adding to a large μ can round a bounded step past maxdx
    if (std::abs(mu_next - mu) > p.maxdx) mu_next = std::nextafter(mu_next, mu);

    finish(mu_next, std::move(status));
    if (std::abs(result.report.filling_next - goal) > std::abs(error_old) + p.root_tol)
      throw std::runtime_error("mu_control_t: the step to μ = " + number(mu_next) + " made the filling error worse");
    return result;
  }

} // namespace nrgljubljana_interface::dmft
