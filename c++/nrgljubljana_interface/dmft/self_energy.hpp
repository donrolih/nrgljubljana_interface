// The self-energy from the correlators of an impurity solve, and the impurity Green's function.

#pragma once

#include "./common.hpp"

namespace nrgljubljana_interface::dmft {

  /**
   * The self-energy from G, F_l, F_r and I: Σ = Σ_H + I − F_l G⁻¹ F_r, per block and mesh point, with the constant
   * Hartree term Σ_H (one matrix per block).
   *
   * F_l, F_r and I must have the structure and mesh of G. Throws std::invalid_argument for a mismatch, a non-finite
   * input, a wrong number or size of Hartree matrices, and a singular G (naming the block and the frequency).
   */
  g_w_t sigma_IFG(g_w_t const &G, g_w_t const &F_l, g_w_t const &F_r, g_w_t const &I, block_matrices_t const &sigma_hartree);

  /// The self-energy Σ = F_l G⁻¹, the estimator without I and F_r. Errors as for sigma_IFG.
  g_w_t sigma_FG(g_w_t const &G, g_w_t const &F_l);

  /**
   * The impurity Green's function from the impurity Dyson equation, G = (G0⁻¹ − Σ)⁻¹ with G0⁻¹ = ω − h_imp − Δ,
   * per block and mesh point, with the constant impurity levels h_imp (one matrix per block).
   *
   * Σ must have the structure and mesh of Δ. Throws std::invalid_argument for a mismatch, a non-finite input, a
   * wrong number or size of h_imp, and a singular ω − h_imp − Δ − Σ (naming the block and the frequency).
   */
  g_w_t impurity_gf(g_w_t const &Delta, g_w_t const &Sigma, block_matrices_t const &h_imp);

} // namespace nrgljubljana_interface::dmft
