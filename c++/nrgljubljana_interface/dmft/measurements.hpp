// Integrals of the spectral matrix A = (G† − G)/(2πi) of a Green's function: occupancy, spectral weight and the
// distance between two spectral functions. Each element is integrated as the Steffen interpolant of its real and
// imaginary parts on the mesh.

#pragma once

#include "./common.hpp"

#include <map>
#include <string>
#include <vector>

namespace nrgljubljana_interface::dmft {

  /// The occupancy of a Green's function
  struct occupancy_t {
    /// N_ij = ∫ f(ω) A_ij(ω) dω per block, Hermitian
    block_matrices_t matrix = {};
    /// The sum over the blocks of degeneracy × Re tr N, e.g. the spin-summed filling
    double total = 0;
  };

  /**
   * The occupancy of G at temperature T (T = 0: the step at ω = 0). degeneracy gives the number of equivalent copies
   * each block stands for (the manifest's [structure] degeneracy), for every block of G and no other.
   *
   * Throws std::invalid_argument for a non-finite G, a negative or non-finite T, and a degeneracy that does not cover
   * exactly the blocks of G with positive counts.
   */
  occupancy_t occupancy(g_w_t const &G, double T, std::map<std::string, long> const &degeneracy, double epsabs = 1e-10,
                        double epsrel = 1e-9);

  /// ∫ A_ij(ω) dω per block, Hermitian: 1 on the diagonal for a normalized G, up to the weight beyond the mesh
  block_matrices_t spectral_weight(g_w_t const &G);

  /// The distance between two spectral functions
  struct distance_t {
    /// Per block and element: ∫ |Re ΔA_ij| dω + ∫ |Im ΔA_ij| dω, which is ∫ |ΔA_ii| dω on the diagonal
    std::vector<nda::matrix<double>> elements = {};
    /// The largest element over all blocks
    double max = 0;
  };

  /// The distance between the spectral functions of a and b, which must have the same structure and mesh
  distance_t spectral_distance(g_w_t const &a, g_w_t const &b);

} // namespace nrgljubljana_interface::dmft
