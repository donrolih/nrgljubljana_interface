// The lattice of the DMFT self-consistency: from the local self-energy to the local lattice Green's function and the
// hybridization for which the impurity reproduces it.

#pragma once

#include "../common.hpp"

#include <complex>
#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace nrgljubljana_interface::dmft {

  /**
   * The lattice protocol. With z = ω + μ − Σ(ω) and a local Σ, the moments are G_n = (1/N) Σ_k h_kⁿ (z − h_k)⁻¹ for
   * n = 0, 1, in the block structure and mesh of Σ. G0 is the local lattice Green's function, and
   * Δ = ω + μ − Σ − G0⁻¹ = G1 G0⁻¹. The on-site energy ⟨h⟩ = (1/N) Σ_k h_k is the high-frequency limit of Δ.
   */
  class lattice_t {
   public:
    virtual ~lattice_t() = default;

    /// (G0, G1)
    virtual std::pair<g_w_t, g_w_t> moments(g_w_t const &Sigma, double mu) const = 0;

    /// G0 alone; by default the first of moments
    virtual g_w_t local_gf(g_w_t const &Sigma, double mu) const;

    /// ⟨h⟩ per block of the given structure
    virtual block_matrices_t onsite_energy(triqs::gfs::gf_struct_t const &structure) const = 0;

    /// The type and the parameters, e.g. to check on resume that the lattice is the same
    [[nodiscard]] virtual std::string description() const = 0;

    /**
     * (G_loc, Δ) with G_loc = G0 and Δ = G1 G0⁻¹, which avoids subtracting two large numbers where G0 is small.
     * Throws std::invalid_argument for a singular G0, naming the block and the frequency.
     */
    std::pair<g_w_t, g_w_t> self_consistency(g_w_t const &Sigma, double mu) const;
  };

  namespace detail {

    /// The scalar moment G_n of one orbital at the mesh points, for σ given at each point
    using scalar_moment_t = std::function<std::vector<std::complex<double>>(std::string const &block, long orbital, std::vector<double> const &mesh,
                                                                            std::vector<std::complex<double>> const &sigma, double mu, unsigned moment)>;

    /**
     * The moments 0 .. max_moment of a lattice given by the scalar transform of each orbital.
     *
     * A block whose off-diagonal Σ is zero is transformed element by element. Otherwise, if shared_orbitals (one DOS
     * for all orbitals of the block), Σ(ω) = U diag(σ_m) U⁻¹ at each point and G_n = U diag(g_n(σ_m)) U⁻¹; an
     * eigenvector matrix with a condition number above 1e8 is an error. Without shared_orbitals an off-diagonal Σ is
     * an error. Σ must be finite.
     */
    std::vector<g_w_t> lattice_moments(g_w_t const &Sigma, double mu, unsigned max_moment, scalar_moment_t const &scalar, bool shared_orbitals,
                                       std::string const &function);

    /// Throws std::invalid_argument unless the floor of Im Σ is finite and at least what NRG Ljubljana accepts
    void check_Sigma_floor(double floor, std::string const &function);

  } // namespace detail

} // namespace nrgljubljana_interface::dmft
