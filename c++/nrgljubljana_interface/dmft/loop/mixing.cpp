#include "./mixing.hpp"

#include <nda/linalg.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace nrgljubljana_interface::dmft {

  namespace detail {

    std::vector<double> integration_weights(std::vector<double> const &points) {
      std::size_t const n = points.size();
      std::vector<double> w(n, 0.0);
      for (std::size_t k = 0; k + 1 < n; ++k) {
        double const half = (points[k + 1] - points[k]) / 2;
        w[k] += half;
        w[k + 1] += half;
      }
      return w;
    }

    double weighted_distance(g_w_t const &a, g_w_t const &b) {
      double sum = 0;
      for (long bl = 0; bl < long(a.size()); ++bl) {
        auto const &x = a[bl].data();
        auto const &y = b[bl].data();
        auto const w  = integration_weights(a[bl].mesh().points());
        for (long k = 0; k < x.extent(0); ++k)
          for (long i = 0; i < x.extent(1); ++i)
            for (long j = 0; j < x.extent(2); ++j) sum += w[k] * std::norm(x(k, i, j) - y(k, i, j));
      }
      return std::sqrt(sum);
    }

    void check_mixer_type(h5::group g, std::string const &name, std::string const &type) {
      auto grp = g.open_group(name);
      std::string stored;
      h5::read(grp, "type", stored);
      if (stored != type) throw std::runtime_error(type + " mixer: the group " + name + " holds a " + stored + " mixer");
    }

  } // namespace detail

  hybridization_mix_t mixer_t::mix_hybridization(long iteration, g_w_t const &Delta, g_w_t const &Delta_out, block_matrices_t const &onsite_energy,
                                                  double Gamma_floor, double mu, double filling, std::optional<double> goal) {
    auto mixed                  = mix({iteration, spectral_part(Delta), spectral_part(Delta_out), mu, filling, goal});
    auto [Delta_mix, projection] = hybridization_from_Gamma(mixed.Gamma_mix, onsite_energy, Gamma_floor);
    return {std::move(Delta_mix), mixed.mu_next, std::move(mixed.report), projection};
  }

  linear_mixer::linear_mixer(double alpha) : alpha_(alpha) {
    if (not(alpha_ > 0 and alpha_ <= 1)) throw std::invalid_argument("linear_mixer: α must satisfy 0 < α ≤ 1, got " + detail::number(alpha_));
  }

  mix_result_t linear_mixer::mix(mix_input_t const &input) {
    detail::check_same_structure(input.Gamma_out, input.Gamma, "linear_mixer", "Gamma_out", "Gamma");
    detail::check_finite(input.Gamma, "linear_mixer", "Gamma");
    detail::check_finite(input.Gamma_out, "linear_mixer", "Gamma_out");

    g_w_t next = input.Gamma;
    for (long b = 0; b < long(next.size()); ++b) {
      auto &out       = next[b].data();
      auto const &raw = input.Gamma_out[b].data();
      auto const &old = input.Gamma[b].data();
      for (long k = 0; k < out.extent(0); ++k)
        for (long i = 0; i < out.extent(1); ++i)
          for (long j = 0; j < out.extent(2); ++j) out(k, i, j) = alpha_ * raw(k, i, j) + (1 - alpha_) * old(k, i, j);
    }
    mix_report_t report;
    report.method        = "linear";
    report.residual_norm = detail::weighted_distance(input.Gamma_out, input.Gamma);
    return {std::move(next), std::nullopt, std::move(report)};
  }

  std::string linear_mixer::description() const {
    std::ostringstream os;
    os.precision(17);
    os << "linear alpha=" << alpha_;
    return os.str();
  }

  void linear_mixer::h5_write(h5::group g, std::string const &name) const {
    auto grp = g.create_group(name);
    h5::write(grp, "type", std::string("linear"));
    h5::write(grp, "alpha", alpha_);
  }

  void linear_mixer::h5_read(h5::group g, std::string const &name) { detail::check_mixer_type(g, name, "linear"); }

  void h5_write(h5::group g, std::string const &name, mix_report_t const &r) {
    auto grp = g.create_group(name);
    h5::write_hdf5_format(grp, r);
    h5::write(grp, "method", r.method);
    h5::write(grp, "residual_norm", r.residual_norm);
    h5::write(grp, "history_size", r.history_size);
    h5::write(grp, "condition", r.condition);
    h5::write(grp, "fallback", long(r.fallback));
  }

  void h5_read(h5::group g, std::string const &name, mix_report_t &r) {
    auto grp = g.open_group(name);
    long fallback = 0;
    h5::read(grp, "method", r.method);
    h5::read(grp, "residual_norm", r.residual_norm);
    h5::read(grp, "history_size", r.history_size);
    h5::read(grp, "condition", r.condition);
    h5::read(grp, "fallback", fallback);
    r.fallback = fallback != 0;
  }

  // ---------------------------------------------------------------- Broyden

  namespace {

    // The scale of the components at each mesh point of a block: √w_k, or 1
    std::vector<double> component_scale(g_w_t const &G, long b, bool weighted) {
      auto const pts = G[b].mesh().points();
      std::vector<double> scale(pts.size(), 1.0);
      if (weighted) {
        auto const w = detail::integration_weights(std::vector<double>(pts.begin(), pts.end()));
        for (std::size_t k = 0; k < w.size(); ++k) scale[k] = std::sqrt(w[k]);
      }
      return scale;
    }

    // The independent real components of a Hermitian Γ: per block and mesh point the diagonal, then Re and Im of the
    // upper triangle
    std::vector<double> pack(g_w_t const &Gamma, bool weighted) {
      std::vector<double> x;
      for (long b = 0; b < long(Gamma.size()); ++b) {
        auto const &d     = Gamma[b].data();
        auto const scale  = component_scale(Gamma, b, weighted);
        long const n      = d.extent(1);
        for (long k = 0; k < d.extent(0); ++k) {
          for (long i = 0; i < n; ++i) x.push_back(d(k, i, i).real() * scale[k]);
          for (long i = 0; i < n; ++i)
            for (long j = i + 1; j < n; ++j) {
              x.push_back(d(k, i, j).real() * scale[k]);
              x.push_back(d(k, i, j).imag() * scale[k]);
            }
        }
      }
      return x;
    }

    g_w_t unpack(std::vector<double> const &x, g_w_t const &like, bool weighted) {
      g_w_t Gamma = like;
      std::size_t pos = 0;
      for (long b = 0; b < long(Gamma.size()); ++b) {
        auto &d          = Gamma[b].data();
        auto const scale = component_scale(Gamma, b, weighted);
        long const n     = d.extent(1);
        for (long k = 0; k < d.extent(0); ++k) {
          for (long i = 0; i < n; ++i) d(k, i, i) = x[pos++] / scale[k];
          for (long i = 0; i < n; ++i)
            for (long j = i + 1; j < n; ++j) {
              std::complex<double> const v{x[pos] / scale[k], x[pos + 1] / scale[k]};
              pos += 2;
              d(k, i, j) = v;
              d(k, j, i) = std::conj(v);
            }
        }
      }
      return Gamma;
    }

    // The smallest eigenvalue of Γ over all blocks and mesh points, and the largest |Γ_ij|
    std::pair<double, double> eigenvalue_range(g_w_t const &Gamma) {
      double lowest = std::numeric_limits<double>::infinity(), peak = 0;
      for (long b = 0; b < long(Gamma.size()); ++b) {
        auto const &d = Gamma[b].data();
        long const n  = d.extent(1);
        for (long k = 0; k < d.extent(0); ++k) {
          nda::matrix<std::complex<double>> m(n, n);
          for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j) {
              m(i, j) = d(k, i, j);
              peak    = std::max(peak, std::abs(m(i, j)));
            }
          lowest = std::min(lowest, n == 1 ? m(0, 0).real() : double(nda::linalg::eigvalsh(m)(0)));
        }
      }
      return {lowest, peak};
    }

    struct step_t {
      std::vector<double> x;
      long pairs       = 0;
      double condition = 1;
    };

    double norm(std::vector<double> const &v) {
      double s = 0;
      for (double a : v) s += a * a;
      return std::sqrt(s);
    }

    // The modified Broyden step from the samples (inputs X, residuals F, chronological) with the factors a
    step_t broyden_step(std::vector<std::vector<double>> const &X, std::vector<std::vector<double>> const &F, std::vector<double> const &a, double w0) {
      std::size_t const N = a.size();
      auto const &x = X.back(), &f = F.back();
      step_t step;
      step.x.resize(N);
      for (std::size_t i = 0; i < N; ++i) step.x[i] = x[i] + a[i] * f[i];

      std::vector<std::vector<double>> delta, correction;
      for (std::size_t s = 0; s + 1 < X.size(); ++s) {
        std::vector<double> dF(N);
        for (std::size_t i = 0; i < N; ++i) dF[i] = F[s + 1][i] - F[s][i];
        double const length = norm(dF);
        double const scale  = std::max({norm(F[s]), norm(F[s + 1]), 1.0});
        if (length <= 10 * std::numeric_limits<double>::epsilon() * scale) continue;
        std::vector<double> u(N), c(N);
        for (std::size_t i = 0; i < N; ++i) {
          u[i] = dF[i] / length;
          c[i] = a[i] * u[i] + (X[s + 1][i] - X[s][i]) / length;
        }
        delta.push_back(std::move(u));
        correction.push_back(std::move(c));
      }
      step.pairs = long(delta.size());
      if (step.pairs == 0) return step;

      long const m = step.pairs;
      nda::matrix<double> system(m, m);
      nda::vector<double> projection(m);
      for (long p = 0; p < m; ++p) {
        double proj = 0;
        for (std::size_t i = 0; i < N; ++i) proj += delta[p][i] * f[i];
        projection(p) = proj;
        for (long q = 0; q <= p; ++q) {
          double g = 0;
          for (std::size_t i = 0; i < N; ++i) g += delta[p][i] * delta[q][i];
          system(p, q) = system(q, p) = g;
        }
        system(p, p) += w0;
      }
      auto const eigenvalues = nda::linalg::eigvalsh(system);
      step.condition         = eigenvalues(m - 1) / eigenvalues(0);
      nda::vector<double> const coefficients = nda::linalg::solve(system, projection);
      for (long p = 0; p < m; ++p)
        for (std::size_t i = 0; i < N; ++i) step.x[i] -= correction[p][i] * coefficients(p);
      return step;
    }

    bool all_finite(std::vector<double> const &v) {
      return std::all_of(v.begin(), v.end(), [](double a) { return std::isfinite(a); });
    }

  } // namespace

  void h5_write(h5::group g, std::string const &name, broyden_params_t const &p) {
    auto grp = g.create_group(name);
    h5::write_hdf5_format(grp, p);
    h5::write(grp, "alpha", p.alpha);
    h5::write(grp, "history", p.history);
    h5::write(grp, "mu", long(p.mu));
    h5::write(grp, "mu_alpha", p.mu_alpha);
    h5::write(grp, "warmup", p.warmup);
    h5::write(grp, "max_mu_step", p.max_mu_step);
    h5::write(grp, "weighted", long(p.weighted));
    h5::write(grp, "w0", p.w0);
    h5::write(grp, "max_condition", p.max_condition);
  }

  void h5_read(h5::group g, std::string const &name, broyden_params_t &p) {
    auto grp      = g.open_group(name);
    long mu       = 0;
    long weighted = 0;
    h5::read(grp, "alpha", p.alpha);
    h5::read(grp, "history", p.history);
    h5::read(grp, "mu", mu);
    h5::read(grp, "mu_alpha", p.mu_alpha);
    h5::read(grp, "warmup", p.warmup);
    h5::read(grp, "max_mu_step", p.max_mu_step);
    h5::read(grp, "weighted", weighted);
    h5::read(grp, "w0", p.w0);
    h5::read(grp, "max_condition", p.max_condition);
    p.mu       = mu != 0;
    p.weighted = weighted != 0;
  }

  broyden_mixer::broyden_mixer(broyden_params_t const &params) : params_(params) {
    auto const &p   = params_;
    std::string pre = "broyden_mixer: ";
    auto positive   = [](double x) { return std::isfinite(x) and x > 0; };
    if (not(p.alpha > 0 and p.alpha <= 1)) throw std::invalid_argument(pre + "alpha must satisfy 0 < alpha ≤ 1, got " + detail::number(p.alpha));
    if (p.history < 1) throw std::invalid_argument(pre + "history must be at least 1, got " + std::to_string(p.history));
    if (not positive(p.mu_alpha)) throw std::invalid_argument(pre + "mu_alpha must be positive, got " + detail::number(p.mu_alpha));
    if (p.warmup < 0) throw std::invalid_argument(pre + "warmup must be >= 0, got " + std::to_string(p.warmup));
    if (not positive(p.max_mu_step)) throw std::invalid_argument(pre + "max_mu_step must be positive, got " + detail::number(p.max_mu_step));
    if (not(std::isfinite(p.w0) and p.w0 >= 0)) throw std::invalid_argument(pre + "w0 must be finite and >= 0, got " + detail::number(p.w0));
    if (not(p.max_condition >= 1)) throw std::invalid_argument(pre + "max_condition must be at least 1, got " + detail::number(p.max_condition));
  }

  mix_result_t broyden_mixer::mix(mix_input_t const &input) {
    auto const &p = params_;
    std::string const pre = "broyden_mixer: ";
    detail::check_same_structure(input.Gamma_out, input.Gamma, "broyden_mixer", "Gamma_out", "Gamma");
    detail::check_finite(input.Gamma, "broyden_mixer", "Gamma");
    detail::check_finite(input.Gamma_out, "broyden_mixer", "Gamma_out");
    if (input.iteration < 1) throw std::invalid_argument(pre + "the iteration must be at least 1");
    if (not history_.empty() and input.iteration != history_.back().iteration + 1)
      throw std::invalid_argument(pre + "iteration " + std::to_string(input.iteration) + " does not follow iteration "
                                  + std::to_string(history_.back().iteration));
    if (p.mu and not input.goal) throw std::invalid_argument(pre + "mixing μ needs the goal of the filling");
    if (p.mu and not(std::isfinite(input.mu) and std::isfinite(input.filling))) throw std::invalid_argument(pre + "μ and the filling must be finite");

    sample_t sample;
    sample.iteration = input.iteration;
    sample.x         = pack(input.Gamma, p.weighted);
    auto const out   = pack(input.Gamma_out, p.weighted);
    sample.F.resize(out.size());
    for (std::size_t i = 0; i < out.size(); ++i) sample.F[i] = out[i] - sample.x[i];
    sample.mu      = input.mu;
    sample.filling = input.filling;
    if (not history_.empty() and history_.back().x.size() != sample.x.size())
      throw std::invalid_argument(pre + "Γ has another structure than in the previous iteration");
    history_.push_back(std::move(sample));
    if (long(history_.size()) > p.history) history_.erase(history_.begin(), history_.end() - p.history);

    bool const with_mu   = moves_mu(input.iteration);
    std::size_t const NG = history_.back().x.size();
    std::vector<std::vector<double>> X, F;
    for (auto const &s : history_) {
      X.push_back(s.x);
      F.push_back(s.F);
      if (with_mu) {
        X.back().push_back(s.mu);
        F.back().push_back(*input.goal - s.filling);
      }
    }
    std::vector<double> a(NG, p.alpha);
    if (with_mu) a.push_back(p.mu_alpha);

    auto linear = [&] {
      std::vector<double> x(a.size());
      for (std::size_t i = 0; i < a.size(); ++i) x[i] = X.back()[i] + a[i] * F.back()[i];
      return x;
    };
    step_t step = broyden_step(X, F, a, p.w0);
    bool fallback = step.pairs > 0 and (not all_finite(step.x) or not(step.condition <= p.max_condition));
    if (fallback) step.x = linear();
    g_w_t Gamma_mix = unpack(std::vector<double>(step.x.begin(), step.x.begin() + long(NG)), input.Gamma, p.weighted);
    if (step.pairs > 0 and not fallback) {
      auto const [lowest, peak] = eigenvalue_range(Gamma_mix);
      if (lowest < -1e-12 * std::max(1.0, peak)) {
        fallback  = true;
        step.x    = linear();
        Gamma_mix = unpack(std::vector<double>(step.x.begin(), step.x.begin() + long(NG)), input.Gamma, p.weighted);
      }
    }

    mix_result_t result;
    result.Gamma_mix = std::move(Gamma_mix);
    if (with_mu) {
      double mu_next = input.mu + std::clamp(step.x.back() - input.mu, -p.max_mu_step, p.max_mu_step);
      if (std::abs(mu_next - input.mu) > p.max_mu_step) mu_next = std::nextafter(mu_next, input.mu);
      result.mu_next = mu_next;
    }
    result.report.method        = fallback ? "linear (fallback)" : "broyden";
    result.report.residual_norm = detail::weighted_distance(input.Gamma_out, input.Gamma);
    result.report.history_size  = long(history_.size());
    result.report.condition     = step.condition;
    result.report.fallback      = fallback;
    return result;
  }

  std::string broyden_mixer::description() const {
    auto const &p = params_;
    std::ostringstream os;
    os.precision(17);
    os << "broyden alpha=" << p.alpha << " history=" << p.history << " mu=" << p.mu << " mu_alpha=" << p.mu_alpha << " warmup=" << p.warmup
       << " max_mu_step=" << p.max_mu_step << " weighted=" << p.weighted << " w0=" << p.w0 << " max_condition=" << p.max_condition;
    return os.str();
  }

  void broyden_mixer::h5_write(h5::group g, std::string const &name) const {
    auto grp = g.create_group(name);
    h5::write(grp, "type", std::string("broyden"));
    nrgljubljana_interface::dmft::h5_write(grp, "params", params_);
    long const n = long(history_.size());
    long const N = n > 0 ? long(history_.front().x.size()) : 0;
    nda::array<double, 2> x(n, N), F(n, N);
    nda::array<long, 1> iterations(n);
    nda::array<double, 1> mu(n), filling(n);
    for (long s = 0; s < n; ++s) {
      iterations(s) = history_[s].iteration;
      mu(s)         = history_[s].mu;
      filling(s)    = history_[s].filling;
      for (long i = 0; i < N; ++i) {
        x(s, i) = history_[s].x[i];
        F(s, i) = history_[s].F[i];
      }
    }
    h5::write(grp, "iterations", iterations);
    h5::write(grp, "x", x);
    h5::write(grp, "F", F);
    h5::write(grp, "mu", mu);
    h5::write(grp, "filling", filling);
  }

  void broyden_mixer::h5_read(h5::group g, std::string const &name) {
    detail::check_mixer_type(g, name, "broyden");
    auto grp = g.open_group(name);
    nda::array<double, 2> x, F;
    nda::array<long, 1> iterations;
    nda::array<double, 1> mu, filling;
    h5::read(grp, "iterations", iterations);
    h5::read(grp, "x", x);
    h5::read(grp, "F", F);
    h5::read(grp, "mu", mu);
    h5::read(grp, "filling", filling);
    history_.clear();
    for (long s = 0; s < iterations.extent(0); ++s) {
      sample_t sample;
      sample.iteration = iterations(s);
      sample.mu        = mu(s);
      sample.filling   = filling(s);
      for (long i = 0; i < x.extent(1); ++i) {
        sample.x.push_back(x(s, i));
        sample.F.push_back(F(s, i));
      }
      history_.push_back(std::move(sample));
    }
    if (long(history_.size()) > params_.history) history_.erase(history_.begin(), history_.end() - params_.history);
  }

} // namespace nrgljubljana_interface::dmft
