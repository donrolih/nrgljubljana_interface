// A lattice given by a tabulated density of states.

#pragma once

#include "./lattice.hpp"

#include <map>
#include <optional>
#include <string>
#include <vector>

namespace nrgljubljana_interface::dmft {

  /// A tabulated density of states: strictly increasing energies and non-negative densities, normalized to 1
  struct dos_t {
    std::vector<double> energy, density;
  };

  /**
   * A lattice given by one tabulated DOS for all orbitals, or one per diagonal element of each block.
   *
   * The moments are Hilbert transforms of the Steffen interpolant of the DOS (NRG Ljubljana numerics), at
   * z = ω + μ − Re Σ + i max(−Im Σ, floor). With one DOS for all orbitals a matrix Σ goes through its eigenvectors;
   * with a DOS per element an off-diagonal Σ is an error. The on-site energy is ∫ ε ρ(ε) dε on the diagonal.
   *
   * Each DOS is checked at construction: at least 3 finite points, strictly increasing energies, ρ ≥ 0, and
   * ∫ ρ dε = 1 within 1e-8. It is not rescaled. Errors are std::invalid_argument.
   */
  class tabulated_lattice_t : public lattice_t {
   public:
    explicit tabulated_lattice_t(dos_t dos, double sigma_floor = 1e-12);
    explicit tabulated_lattice_t(std::map<std::string, std::vector<dos_t>> per_element, double sigma_floor = 1e-12);

    std::pair<g_w_t, g_w_t> moments(g_w_t const &Sigma, double mu) const override;
    g_w_t local_gf(g_w_t const &Sigma, double mu) const override;
    block_matrices_t onsite_energy(triqs::gfs::gf_struct_t const &structure) const override;
    /// The kind (shared or per element), the number of points and a checksum of each table, and the floor
    [[nodiscard]] std::string description() const override;

    [[nodiscard]] double sigma_floor() const { return floor_; }

    /// A DOS from a file of two columns, energy and density; blank lines and lines starting with # are skipped
    static dos_t read_dos(std::string const &path);

   private:
    std::optional<dos_t> shared_;
    std::map<std::string, std::vector<dos_t>> per_element_;
    std::map<std::string, std::vector<double>> centroid_; // per element; key "" for the shared DOS
    double floor_;

    [[nodiscard]] dos_t const &dos(std::string const &block, long orbital) const;
    void check_structure(std::vector<std::string> const &names, std::vector<long> const &sizes, std::string const &function) const;
    [[nodiscard]] std::vector<g_w_t> transform(g_w_t const &Sigma, double mu, unsigned max_moment, std::string const &function) const;
  };

} // namespace nrgljubljana_interface::dmft
