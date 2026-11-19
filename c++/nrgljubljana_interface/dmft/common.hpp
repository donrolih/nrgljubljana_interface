// What the DMFT functions share: one matrix per block, and the checks of their inputs.

#pragma once

#include "../types.hpp"

#include <nda/nda.hpp>

#include <complex>
#include <string>
#include <vector>

namespace nrgljubljana_interface::dmft {

  /// One matrix per block, in the order of the blocks
  using block_matrices_t = std::vector<nda::matrix<std::complex<double>>>;

  namespace detail {

    /// A number for an error message
    std::string number(double x);

    /// Throws std::invalid_argument "<function>: <what> is not finite in block B, element (i, j), at ω = x"
    void check_finite(g_w_t const &F, std::string const &function, std::string const &what = "the function");

    /// Throws std::invalid_argument unless G has the block names, block sizes and meshes of reference
    void check_same_structure(g_w_t const &G, g_w_t const &reference, std::string const &function, std::string const &what,
                              std::string const &reference_name);

    /// Throws std::invalid_argument unless there is one finite matrix per block of F, of the size of the block
    void check_block_matrices(g_w_t const &F, block_matrices_t const &matrices, std::string const &function, std::string const &what);

    /// The Green's function with the structure and mesh of F, equal to matrices[b] at every point of block b
    g_w_t constant_gf(g_w_t const &F, block_matrices_t const &matrices);

  } // namespace detail

} // namespace nrgljubljana_interface::dmft
