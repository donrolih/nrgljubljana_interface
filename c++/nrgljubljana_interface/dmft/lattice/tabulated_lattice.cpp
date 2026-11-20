#include "./tabulated_lattice.hpp"

#include "../../nrg/integrate.hpp"
#include "../../nrg/lattice.hpp"
#include "../../nrg/output.hpp"

#include <nda/linalg.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace nrgljubljana_interface::dmft {

  namespace {

    // A DOS is normalized if its integral is within this of 1
    constexpr double normalization_tolerance = 1e-8;

    nrg::dos_table_t table(dos_t const &d) { return {d.energy, d.density}; }

    void check_dos(dos_t const &d, std::string const &where, std::string const &function = "tabulated_lattice_t", bool normalized = true) {
      auto const prefix = function + ": the DOS" + where;
      if (d.energy.size() != d.density.size())
        throw std::invalid_argument(prefix + " has " + std::to_string(d.energy.size()) + " energies and " + std::to_string(d.density.size()) + " densities");
      if (d.energy.size() < 3) throw std::invalid_argument(prefix + " has fewer than 3 points");
      for (std::size_t k = 0; k < d.energy.size(); ++k) {
        if (not std::isfinite(d.energy[k]) or not std::isfinite(d.density[k])) throw std::invalid_argument(prefix + " is not finite at point " + std::to_string(k));
        if (d.density[k] < 0)
          throw std::invalid_argument(prefix + " is negative at ε = " + detail::number(d.energy[k]) + " (" + detail::number(d.density[k]) + ")");
        if (k > 0 and not(d.energy[k] > d.energy[k - 1]))
          throw std::invalid_argument(prefix + ": the energies are not strictly increasing at point " + std::to_string(k));
      }
      if (not normalized) return;
      double const norm = nrg::integrate(d.energy, d.density);
      if (not(std::abs(norm - 1) <= normalization_tolerance))
        throw std::invalid_argument(prefix + " is not normalized: its integral is " + detail::number(norm) + ", expected 1 within 1e-8");
    }

    double centroid(dos_t const &d) { return nrg::integrate_moment(table(d), 1); }

    // FNV-1a over the bytes of the energies and densities: a table changed in any value changes it
    std::uint64_t checksum(dos_t const &d) {
      std::uint64_t h = 1469598103934665603ULL;
      auto add        = [&](std::vector<double> const &v) {
        auto const *p = reinterpret_cast<unsigned char const *>(v.data());
        for (std::size_t i = 0; i < v.size() * sizeof(double); ++i) h = (h ^ p[i]) * 1099511628211ULL;
      };
      add(d.energy);
      add(d.density);
      return h;
    }

    std::string describe(dos_t const &d) {
      std::ostringstream os;
      os << d.energy.size() << " points, checksum " << std::hex << checksum(d);
      return os.str();
    }

  } // namespace

  tabulated_lattice_t::tabulated_lattice_t(dos_t dos, double Sigma_floor) : shared_(std::move(dos)), floor_(Sigma_floor) {
    detail::check_Sigma_floor(floor_, "tabulated_lattice_t");
    check_dos(*shared_, "");
    centroid_[""] = {centroid(*shared_)};
  }

  tabulated_lattice_t::tabulated_lattice_t(std::map<std::string, std::vector<dos_t>> per_element, double Sigma_floor)
     : per_element_(std::move(per_element)), floor_(Sigma_floor) {
    detail::check_Sigma_floor(floor_, "tabulated_lattice_t");
    if (per_element_.empty()) throw std::invalid_argument("tabulated_lattice_t: no blocks");
    for (auto const &[name, list] : per_element_) {
      if (list.empty()) throw std::invalid_argument("tabulated_lattice_t: no DOS for block " + name);
      for (std::size_t i = 0; i < list.size(); ++i) {
        check_dos(list[i], " of block " + name + ", orbital " + std::to_string(i));
        centroid_[name].push_back(centroid(list[i]));
      }
    }
  }

  std::string tabulated_lattice_t::description() const {
    std::ostringstream os;
    os.precision(17);
    if (shared_)
      os << "tabulated shared: " << describe(*shared_);
    else {
      os << "tabulated per element:";
      for (auto const &[name, list] : per_element_)
        for (std::size_t i = 0; i < list.size(); ++i) os << " " << name << "[" << i << "]: " << describe(list[i]) << ";";
    }
    os << " Sigma_floor=" << floor_;
    return os.str();
  }

  dos_t const &tabulated_lattice_t::dos(std::string const &block, long orbital) const { return shared_ ? *shared_ : per_element_.at(block).at(orbital); }

  void tabulated_lattice_t::check_structure(std::vector<std::string> const &names, std::vector<long> const &sizes, std::string const &function) const {
    if (shared_) return;
    for (std::size_t b = 0; b < names.size(); ++b) {
      auto const it = per_element_.find(names[b]);
      if (it == per_element_.end()) throw std::invalid_argument(function + ": no DOS for block " + names[b]);
      if (long(it->second.size()) != sizes[b])
        throw std::invalid_argument(function + ": block " + names[b] + " has " + std::to_string(sizes[b]) + " orbitals, the DOS map "
                                    + std::to_string(it->second.size()));
    }
    for (auto const &[name, list] : per_element_)
      if (std::find(names.begin(), names.end(), name) == names.end()) throw std::invalid_argument(function + ": a DOS for block " + name + ", not in the structure");
  }

  std::vector<g_w_t> tabulated_lattice_t::transform(g_w_t const &Sigma, double mu, unsigned max_moment, std::string const &function) const {
    std::vector<long> sizes;
    for (long b = 0; b < long(Sigma.size()); ++b) sizes.push_back(Sigma[b].data().extent(1));
    check_structure(Sigma.block_names(), sizes, function);
    auto const floor = floor_;
    auto scalar = [this, floor](std::string const &block, long orbital, std::vector<double> const &mesh, std::vector<std::complex<double>> const &sigma,
                                double m, unsigned moment) {
      return nrg::hilbert_moment(table(dos(block, orbital)), mesh, sigma, m, floor, moment, detail::workers());
    };
    return detail::lattice_moments(Sigma, mu, max_moment, scalar, shared_.has_value(), function);
  }

  std::pair<g_w_t, g_w_t> tabulated_lattice_t::moments(g_w_t const &Sigma, double mu) const {
    auto m = transform(Sigma, mu, 1, "tabulated_lattice_t");
    return {std::move(m[0]), std::move(m[1])};
  }

  g_w_t tabulated_lattice_t::local_gf(g_w_t const &Sigma, double mu) const { return std::move(transform(Sigma, mu, 0, "tabulated_lattice_t")[0]); }

  block_matrices_t tabulated_lattice_t::onsite_energy(triqs::gfs::gf_struct_t const &structure) const {
    std::vector<std::string> names;
    std::vector<long> sizes;
    for (auto const &[name, n] : structure) {
      names.push_back(name);
      sizes.push_back(n);
    }
    check_structure(names, sizes, "tabulated_lattice_t");
    block_matrices_t result;
    for (std::size_t b = 0; b < names.size(); ++b) {
      nda::matrix<std::complex<double>> m(sizes[b], sizes[b]);
      m = 0;
      for (long i = 0; i < sizes[b]; ++i) m(i, i) = shared_ ? centroid_.at("")[0] : centroid_.at(names[b])[i];
      result.push_back(std::move(m));
    }
    return result;
  }

  dos_t tabulated_lattice_t::read_dos(std::string const &path) {
    std::istringstream in(nrg::read_text_file(path));
    dos_t d;
    long number = 0;
    for (std::string line; std::getline(in, line);) {
      ++number;
      auto const first = line.find_first_not_of(" \t\r");
      if (first == std::string::npos or line[first] == '#') continue;
      std::istringstream fields(line);
      double e = 0, rho = 0;
      std::string extra;
      if (not(fields >> e >> rho) or (fields >> extra))
        throw std::invalid_argument("read_dos: " + path + ", line " + std::to_string(number) + ": expected two numbers, energy and density");
      d.energy.push_back(e);
      d.density.push_back(rho);
    }
    return d;
  }

  namespace {

    using dcomplex = std::complex<double>;

    // Largest condition number of the eigenvector matrix of z accepted by the matrix transform
    constexpr double max_condition = 1e8;

    // G_n(z) for each z, through NRG Ljubljana with ω = Re z, μ = 0 and Σ = −i |Im z|, and G_n(z̄) = conj G_n(z). The
    // frequencies of one call must increase strictly: the points are sorted by Re z, and the k-th of the points
    // with the same Re z goes to the k-th call.
    std::vector<dcomplex> transform_points(dos_t const &dos, std::vector<dcomplex> const &z, long moment) {
      if (moment < 0) throw std::invalid_argument("hilbert_transform: the moment must be non-negative, got " + std::to_string(moment));
      check_dos(dos, "", "hilbert_transform", false);
      for (std::size_t k = 0; k < z.size(); ++k)
        if (not std::isfinite(z[k].real()) or not std::isfinite(z[k].imag()))
          throw std::invalid_argument("hilbert_transform: z is not finite at index " + std::to_string(k));

      std::vector<std::size_t> order(z.size());
      for (std::size_t k = 0; k < z.size(); ++k) order[k] = k;
      std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return z[a].real() < z[b].real(); });
      std::vector<std::vector<std::size_t>> passes;
      for (std::size_t i = 0, repeat = 0; i < order.size(); ++i) {
        repeat = (i > 0 and z[order[i]].real() == z[order[i - 1]].real()) ? repeat + 1 : 0;
        if (passes.size() <= repeat) passes.emplace_back();
        passes[repeat].push_back(order[i]);
      }

      double const floor = std::sqrt(std::numeric_limits<double>::min());
      std::vector<dcomplex> g(z.size());
      for (auto const &pass : passes) {
        std::vector<double> mesh;
        std::vector<dcomplex> sigma;
        for (auto k : pass) {
          mesh.push_back(z[k].real());
          sigma.emplace_back(0, -std::abs(z[k].imag()));
        }
        auto const values = nrg::hilbert_moment(table(dos), mesh, sigma, 0, floor, unsigned(moment), detail::workers());
        for (std::size_t i = 0; i < pass.size(); ++i) g[pass[i]] = z[pass[i]].imag() < 0 ? std::conj(values[i]) : values[i];
      }
      return g;
    }

    double frobenius(nda::matrix<dcomplex> const &m) {
      double s = 0;
      for (long i = 0; i < m.extent(0); ++i)
        for (long j = 0; j < m.extent(1); ++j) s += std::norm(m(i, j));
      return std::sqrt(s);
    }

  } // namespace

  dcomplex hilbert_transform(dos_t const &dos, dcomplex z, long moment) { return transform_points(dos, {z}, moment)[0]; }

  nda::array<dcomplex, 1> hilbert_transform(dos_t const &dos, nda::array<dcomplex, 1> const &z, long moment) {
    std::vector<dcomplex> const points(z.begin(), z.end());
    auto const g = transform_points(dos, points, moment);
    nda::array<dcomplex, 1> result(z.extent(0));
    for (long k = 0; k < z.extent(0); ++k) result(k) = g[k];
    return result;
  }

  nda::matrix<dcomplex> hilbert_transform(dos_t const &dos, nda::matrix<dcomplex> const &z, long moment) {
    long const n = z.extent(0);
    if (z.extent(1) != n)
      throw std::invalid_argument("hilbert_transform: z must be a square matrix, got " + std::to_string(n) + "x" + std::to_string(z.extent(1)));
    if (n == 0) return nda::matrix<dcomplex>(0, 0);
    for (long i = 0; i < n; ++i)
      for (long j = 0; j < n; ++j)
        if (not std::isfinite(z(i, j).real()) or not std::isfinite(z(i, j).imag()))
          throw std::invalid_argument("hilbert_transform: z is not finite at (" + std::to_string(i) + ", " + std::to_string(j) + ")");
    auto [lambda, V] = nda::linalg::eig(z);
    nda::matrix<dcomplex> U = V;
    if (nda::linalg::det(U) == dcomplex{0}) throw std::invalid_argument("hilbert_transform: z is not diagonalizable");
    nda::matrix<dcomplex> const U_inv = nda::linalg::inv(U);
    double const condition            = frobenius(U) * frobenius(U_inv);
    if (not(condition <= max_condition))
      throw std::invalid_argument("hilbert_transform: the eigenvectors of z are ill-conditioned (condition number " + detail::number(condition) + ")");
    std::vector<dcomplex> const eigenvalues(lambda.begin(), lambda.end());
    auto const g = transform_points(dos, eigenvalues, moment);
    nda::matrix<dcomplex> diag(n, n);
    diag = 0;
    for (long m = 0; m < n; ++m) diag(m, m) = g[m];
    return U * diag * U_inv;
  }

} // namespace nrgljubljana_interface::dmft
