#include "./convergence.hpp"

#include "../measurements.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace nrgljubljana_interface::dmft {

  namespace {

    using detail::number;

    bool positive(double x) { return std::isfinite(x) and x > 0; }

  } // namespace

  void h5_write(h5::group g, std::string const &name, convergence_params_t const &p) {
    auto grp = g.create_group(name);
    h5::write_hdf5_format(grp, p);
    h5::write(grp, "eps", p.eps);
    h5::write(grp, "window", p.window);
    h5::write(grp, "min_iter", p.min_iter);
    h5::write(grp, "max_iter", p.max_iter);
    h5::write(grp, "filling_tol", p.filling_tol);
    h5::write(grp, "mu_tol", p.mu_tol);
    if (p.eps_imp_loc) h5::write(grp, "eps_imp_loc", *p.eps_imp_loc);
  }

  void h5_read(h5::group g, std::string const &name, convergence_params_t &p) {
    auto grp = g.open_group(name);
    h5::read(grp, "eps", p.eps);
    h5::read(grp, "window", p.window);
    h5::read(grp, "min_iter", p.min_iter);
    h5::read(grp, "max_iter", p.max_iter);
    h5::read(grp, "filling_tol", p.filling_tol);
    h5::read(grp, "mu_tol", p.mu_tol);
    p.eps_imp_loc = {};
    if (grp.has_key("eps_imp_loc")) {
      double e = 0;
      h5::read(grp, "eps_imp_loc", e);
      p.eps_imp_loc = e;
    }
  }

  convergence_t::convergence_t(convergence_params_t const &params) : params_(params) {
    auto const &p          = params_;
    std::string const pre = "convergence_t: ";
    if (not positive(p.eps)) throw std::invalid_argument(pre + "eps must be positive, got " + number(p.eps));
    if (p.window < 1) throw std::invalid_argument(pre + "window must be at least 1, got " + std::to_string(p.window));
    if (p.min_iter < 1 or p.max_iter < p.min_iter)
      throw std::invalid_argument(pre + "need 1 ≤ min_iter ≤ max_iter, got " + std::to_string(p.min_iter) + " and " + std::to_string(p.max_iter));
    if (not positive(p.filling_tol)) throw std::invalid_argument(pre + "filling_tol must be positive, got " + number(p.filling_tol));
    if (not positive(p.mu_tol)) throw std::invalid_argument(pre + "mu_tol must be positive, got " + number(p.mu_tol));
    if (p.eps_imp_loc and not positive(*p.eps_imp_loc)) throw std::invalid_argument(pre + "eps_imp_loc must be positive, got " + number(*p.eps_imp_loc));
  }

  convergence_result_t convergence_t::decide(convergence_input_t const &in) {
    auto const &p                = params_;
    std::string const function = "convergence_t";
    if (in.iteration < 1) throw std::invalid_argument(function + ": the iteration must be at least 1, got " + std::to_string(in.iteration));
    if (last_iteration_ != 0 and in.iteration != last_iteration_ + 1)
      throw std::invalid_argument(function + ": iteration " + std::to_string(in.iteration) + " does not follow iteration " + std::to_string(last_iteration_));
    if (in.iteration > 1 and not in.G_loc_previous)
      throw std::invalid_argument(function + ": iteration " + std::to_string(in.iteration) + " needs the G_loc of the previous iteration");
    detail::check_same_structure(in.G_imp, in.G_loc, function, "G_imp", "G_loc");
    if (not std::isfinite(in.filling) or not std::isfinite(in.mu_step)) throw std::invalid_argument(function + ": the filling and the μ step must be finite");

    convergence_result_t r;
    if (in.G_loc_previous) {
      detail::check_same_structure(*in.G_loc_previous, in.G_loc, function, "G_loc_previous", "G_loc");
      r.distance = spectral_distance(in.G_loc, *in.G_loc_previous).max;
      history_.emplace_back(in.iteration, *r.distance);
    }
    last_iteration_    = in.iteration;
    r.distance_imp_loc = spectral_distance(in.G_imp, in.G_loc).max;

    long const available = std::min<long>(long(history_.size()), p.window);
    for (long i = long(history_.size()) - available; i < long(history_.size()); ++i) r.window_max = std::max(r.window_max, history_[i].second);
    r.window_samples = available;
    r.window_ok      = available == p.window and r.window_max < p.eps;
    r.iterations_ok  = in.iteration >= p.min_iter and in.iteration <= p.max_iter;
    r.filling_ok     = not in.goal or std::abs(in.filling - *in.goal) <= p.filling_tol;
    r.mu_ok          = not in.goal or std::abs(in.mu_step) <= p.mu_tol;
    r.imp_loc_ok     = not p.eps_imp_loc or r.distance_imp_loc <= *p.eps_imp_loc;

    if (r.window_ok and r.iterations_ok and r.filling_ok and r.mu_ok and r.imp_loc_ok)
      r.decision = "converged";
    else if (in.iteration >= p.max_iter)
      r.decision = "stop";
    else
      r.decision = "continue";
    return r;
  }

  void convergence_t::h5_write(h5::group g, std::string const &name) const {
    auto grp = g.create_group(name);
    dmft::h5_write(grp, "params", params_);
    std::vector<long> iterations;
    std::vector<double> distances;
    for (auto const &[it, d] : history_) {
      iterations.push_back(it);
      distances.push_back(d);
    }
    h5::write(grp, "iterations", iterations);
    h5::write(grp, "distances", distances);
    h5::write(grp, "last_iteration", last_iteration_);
  }

  void convergence_t::h5_read(h5::group g, std::string const &name) {
    auto grp = g.open_group(name);
    std::vector<long> iterations;
    std::vector<double> distances;
    h5::read(grp, "iterations", iterations);
    h5::read(grp, "distances", distances);
    if (iterations.size() != distances.size()) throw std::runtime_error("convergence_t: the stored history is corrupt");
    h5::read(grp, "last_iteration", last_iteration_);
    history_.clear();
    for (std::size_t i = 0; i < iterations.size(); ++i) history_.emplace_back(iterations[i], distances[i]);
  }

  void h5_write(h5::group g, std::string const &name, convergence_result_t const &r) {
    auto grp = g.create_group(name);
    h5::write_hdf5_format(grp, r);
    h5::write(grp, "decision", r.decision);
    if (r.distance) h5::write(grp, "distance", *r.distance);
    h5::write(grp, "window_max", r.window_max);
    h5::write(grp, "window_samples", r.window_samples);
    h5::write(grp, "distance_imp_loc", r.distance_imp_loc);
    h5::write(grp, "window_ok", long(r.window_ok));
    h5::write(grp, "iterations_ok", long(r.iterations_ok));
    h5::write(grp, "filling_ok", long(r.filling_ok));
    h5::write(grp, "mu_ok", long(r.mu_ok));
    h5::write(grp, "imp_loc_ok", long(r.imp_loc_ok));
  }

  void h5_read(h5::group g, std::string const &name, convergence_result_t &r) {
    auto grp      = g.open_group(name);
    auto read_bool = [&](char const *key, bool &b) {
      long v = 0;
      h5::read(grp, key, v);
      b = v != 0;
    };
    h5::read(grp, "decision", r.decision);
    r.distance = {};
    if (grp.has_key("distance")) {
      double d = 0;
      h5::read(grp, "distance", d);
      r.distance = d;
    }
    h5::read(grp, "window_max", r.window_max);
    h5::read(grp, "window_samples", r.window_samples);
    h5::read(grp, "distance_imp_loc", r.distance_imp_loc);
    read_bool("window_ok", r.window_ok);
    read_bool("iterations_ok", r.iterations_ok);
    read_bool("filling_ok", r.filling_ok);
    read_bool("mu_ok", r.mu_ok);
    read_bool("imp_loc_ok", r.imp_loc_ok);
  }

} // namespace nrgljubljana_interface::dmft
