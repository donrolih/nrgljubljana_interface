// The Bethe lattice: a semicircular density of states, with the moments in closed form.

#pragma once

#include "./lattice.hpp"

namespace nrgljubljana_interface::dmft {

  /**
   * The Bethe lattice of half-bandwidth D, ρ(ε) = 2 √(D² − ε²) / (π D²), the same for every orbital.
   *
   * G0(ζ) = 2 / (ζ + √(ζ − D) √(ζ + D)) with principal roots, the physical branch for Im ζ > 0; G1 = (D²/4) G0²; the
   * on-site energy is 0. ζ = ω + μ − Re σ + i max(−Im σ, floor), as in the tabulated lattice; a matrix Σ goes through its
   * eigenvectors.
   */
  class bethe_lattice_t : public lattice_t {
   public:
    /// Throws std::invalid_argument for a half-bandwidth that is not positive and finite, and for an invalid floor
    explicit bethe_lattice_t(double half_bandwidth, double Sigma_floor = 1e-12);

    std::pair<g_w_t, g_w_t> moments(g_w_t const &Sigma, double mu) const override;
    g_w_t local_gf(g_w_t const &Sigma, double mu) const override;
    block_matrices_t onsite_energy(triqs::gfs::gf_struct_t const &structure) const override;
    [[nodiscard]] std::string description() const override;

    [[nodiscard]] double half_bandwidth() const { return D_; }
    [[nodiscard]] double Sigma_floor() const { return floor_; }

   private:
    double D_, floor_;
    [[nodiscard]] detail::scalar_moment_t scalar() const;
  };

} // namespace nrgljubljana_interface::dmft
