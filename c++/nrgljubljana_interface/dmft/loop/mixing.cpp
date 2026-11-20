#include "./mixing.hpp"

#include <cmath>
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

} // namespace nrgljubljana_interface::dmft
