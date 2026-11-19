// Causality of retarded functions F(ω) = F_∞ + ∫ ρ(x)/(ω − x + i0) dx with ρ(x) ⪰ 0, per block: Δ, Σ and G.
//
// H = (F + F†)/2 is the Hermitian part and S = (F† − F)/(2i) the spectral part, so F = H − iS. F is causal when
// S(ω) ⪰ 0 at every mesh point and H = F_∞ + KK(−S), with the Kramers-Kronig transform KK applied to the real and
// imaginary parts of each element. F_∞ is ⟨h⟩ for Δ, the Hartree term Σ_H for Σ and 0 for G.

#pragma once

#include "./common.hpp"

#include <utility>

namespace nrgljubljana_interface::dmft {

  /// The spectral part S = (F† − F)/(2i) of each block: Γ for Δ, −Im Σ for a 1x1 Σ, π A for G. Exactly Hermitian.
  g_w_t spectral_part(g_w_t const &F);

  /// What causal_projection changed
  struct causality_report_t {
    /// Number of (block, mesh point) pairs at which an eigenvalue of S was raised
    long projected_points = 0;
    /// The largest amount by which an eigenvalue was raised (floor − λ); 0 if none
    double largest_violation = 0;

    bool operator==(causality_report_t const &) const = default;
  };

  /**
   * Raise the eigenvalues of the spectral part below floor to it, keeping the eigenvectors (S' = S + δS), and
   * correct the Hermitian part by the Kramers-Kronig transform of the change: F' = H + KK(−δS) − iS'.
   *
   * For an analytic F the result is the analytic function with spectral part S' and the same F_∞, which is not
   * needed. An eigenvalue within 1e-14 times the largest |F_ij(ω)| of the floor counts as on it, so a second
   * projection changes nothing. A block in which nothing is raised is returned unchanged.
   *
   * Throws std::invalid_argument for a floor that is negative or not finite, and for a non-finite F (naming the
   * block, the element and the frequency).
   */
  std::pair<g_w_t, causality_report_t> causal_projection(g_w_t const &F, double floor);

  /**
   * The Hermitian part rebuilt from the spectral part and the high-frequency limit: F' = F_∞ + KK(−S) − iS.
   *
   * For a function whose Hermitian part is missing or not trusted, e.g. Δ after mixing Γ. F_∞ has one matrix per
   * block of the size of the block, Hermitian to 1e-12 relative; it is symmetrized. Throws std::invalid_argument
   * otherwise, and for a non-finite F.
   */
  g_w_t rebuild_hermitian_part(g_w_t const &F, block_matrices_t const &high_frequency_limit);

  /**
   * Estimate of the high-frequency limit F_∞: H − KK(−S) averaged over the two ends of the mesh, symmetrized.
   *
   * Exact up to the accuracy of KK when ρ vanishes at and beyond the ends of the mesh. Throws std::invalid_argument
   * for a non-finite F.
   */
  block_matrices_t high_frequency_limit(g_w_t const &F);

} // namespace nrgljubljana_interface::dmft
