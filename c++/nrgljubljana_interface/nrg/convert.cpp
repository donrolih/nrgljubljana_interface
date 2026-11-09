#include "./convert.hpp"

#include <complex>
#include <cstddef>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace nrgljubljana_interface::nrg {

  namespace {

    double get(std::complex<double> z, part_t part) { return part == part_t::real ? z.real() : z.imag(); }

    void set(std::complex<double> &z, double v, part_t part) {
      if (part == part_t::real)
        z.real(v);
      else
        z.imag(v);
    }

    void check_table(NRG::Numerics::Table const &t, triqs::mesh::refreq_log const &mesh) {
      if (t.y.size() != t.x.size())
        throw std::invalid_argument("Table has " + std::to_string(t.x.size()) + " points but " + std::to_string(t.y.size()) + " values");
      check_on_mesh(t.x, mesh);
    }

    void check_element(long i, long j, long n_rows, long n_cols) {
      if (i < 0 or i >= n_rows or j < 0 or j >= n_cols) {
        std::ostringstream os;
        os << "Element (" << i << ", " << j << ") is outside the " << n_rows << "x" << n_cols << " target";
        throw std::out_of_range(os.str());
      }
    }

  } // namespace

  std::vector<double> mesh_points(triqs::mesh::refreq_log const &mesh) { return mesh.points(); }

  void check_on_mesh(std::vector<double> const &x, triqs::mesh::refreq_log const &mesh) {
    auto const &pts = mesh.points();
    if (x.size() != pts.size())
      throw std::invalid_argument("Table has " + std::to_string(x.size()) + " points, the mesh " + std::to_string(pts.size()));
    for (std::size_t k = 0; k < x.size(); ++k) {
      if (x[k] != pts[k]) {
        std::ostringstream os;
        os.precision(std::numeric_limits<double>::max_digits10);
        os << "Table point " << k << " is " << x[k] << ", the mesh point " << pts[k];
        throw std::invalid_argument(os.str());
      }
    }
  }

  NRG::Numerics::Table to_table(triqs::gfs::gf_const_view<triqs::mesh::refreq_log, triqs::gfs::matrix_valued> g, long i, long j,
                                part_t part) {
    auto const &d = g.data();
    check_element(i, j, d.extent(1), d.extent(2));
    NRG::Numerics::Table t{mesh_points(g.mesh()), {}};
    t.y.reserve(t.x.size());
    for (long k = 0; k < d.extent(0); ++k) t.y.push_back(get(d(k, i, j), part));
    return t;
  }

  NRG::Numerics::Table to_table(triqs::gfs::gf_const_view<triqs::mesh::refreq_log, triqs::gfs::scalar_valued> g, part_t part) {
    auto const &d = g.data();
    NRG::Numerics::Table t{mesh_points(g.mesh()), {}};
    t.y.reserve(t.x.size());
    for (long k = 0; k < d.extent(0); ++k) t.y.push_back(get(d(k), part));
    return t;
  }

  void from_table(NRG::Numerics::Table const &t, triqs::gfs::gf_view<triqs::mesh::refreq_log, triqs::gfs::matrix_valued> g, long i, long j,
                  part_t part) {
    auto d = g.data();
    check_element(i, j, d.extent(1), d.extent(2));
    check_table(t, g.mesh());
    for (long k = 0; k < d.extent(0); ++k) set(d(k, i, j), t.y[k], part);
  }

  void from_table(NRG::Numerics::Table const &t, triqs::gfs::gf_view<triqs::mesh::refreq_log, triqs::gfs::scalar_valued> g, part_t part) {
    auto d = g.data();
    check_table(t, g.mesh());
    for (long k = 0; k < d.extent(0); ++k) set(d(k), t.y[k], part);
  }

  NRG::Numerics::SpectrumBins to_nrg(spectrum_bins_t const &s, part_t part) {
    if (part == part_t::imag and not s.is_complex) throw std::invalid_argument("The imaginary part was requested from a real spectrum");
    auto const &w = part == part_t::real ? s.weight : s.weight_imag;
    if (w.size() != s.frequency.size())
      throw std::invalid_argument("Spectrum has " + std::to_string(s.frequency.size()) + " frequencies but " + std::to_string(w.size())
                                  + " weights");
    NRG::Numerics::SpectrumBins bins;
    bins.reserve(w.size());
    for (std::size_t k = 0; k < w.size(); ++k) bins.push_back({s.frequency[k], w[k]});
    return bins;
  }

} // namespace nrgljubljana_interface::nrg
