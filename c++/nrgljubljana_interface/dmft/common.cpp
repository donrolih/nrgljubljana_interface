#include "./common.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace nrgljubljana_interface::dmft::detail {

  namespace {

    bool finite(std::complex<double> z) { return std::isfinite(z.real()) and std::isfinite(z.imag()); }

    std::string size(long n, long m) { return std::to_string(n) + "x" + std::to_string(m); }

  } // namespace

  std::string number(double x) {
    std::ostringstream os;
    os.precision(12);
    os << x;
    return os.str();
  }

  std::size_t workers() { return static_cast<std::size_t>(std::max(1, omp_get_max_threads())); }

  void check_finite(g_w_t const &F, std::string const &function, std::string const &what) {
    for (long b = 0; b < long(F.size()); ++b) {
      auto const &d   = F[b].data();
      auto const &pts = F[b].mesh().points();
      for (long k = 0; k < d.extent(0); ++k)
        for (long i = 0; i < d.extent(1); ++i)
          for (long j = 0; j < d.extent(2); ++j)
            if (not finite(d(k, i, j)))
              throw std::invalid_argument(function + ": " + what + " is not finite in block " + F.block_names()[b] + ", element ("
                                          + std::to_string(i) + ", " + std::to_string(j) + "), at ω = " + number(pts[k]));
    }
  }

  void check_same_structure(g_w_t const &G, g_w_t const &reference, std::string const &function, std::string const &what,
                            std::string const &reference_name) {
    auto const prefix = function + ": " + what + " ";
    if (G.block_names() != reference.block_names()) throw std::invalid_argument(prefix + "does not have the blocks of " + reference_name);
    for (long b = 0; b < long(G.size()); ++b) {
      auto const &name = G.block_names()[b];
      auto const &d    = G[b].data();
      auto const &r    = reference[b].data();
      if (d.extent(1) != r.extent(1) or d.extent(2) != r.extent(2))
        throw std::invalid_argument(prefix + "block " + name + " is " + size(d.extent(1), d.extent(2)) + ", in " + reference_name + " "
                                    + size(r.extent(1), r.extent(2)));
      if (not(G[b].mesh() == reference[b].mesh())) throw std::invalid_argument(prefix + "block " + name + " is not on the mesh of " + reference_name);
    }
  }

  void check_block_matrices(g_w_t const &F, block_matrices_t const &matrices, std::string const &function, std::string const &what) {
    auto const prefix = function + ": " + what;
    if (long(matrices.size()) != long(F.size()))
      throw std::invalid_argument(prefix + " has " + std::to_string(matrices.size()) + " matrices, the function " + std::to_string(F.size())
                                  + " blocks");
    for (long b = 0; b < long(F.size()); ++b) {
      auto const &M    = matrices[b];
      auto const &name = F.block_names()[b];
      long const n     = F[b].data().extent(1);
      if (M.extent(0) != n or M.extent(1) != n)
        throw std::invalid_argument(prefix + " of block " + name + " is " + size(M.extent(0), M.extent(1)) + ", the block " + size(n, n));
      for (long i = 0; i < n; ++i)
        for (long j = 0; j < n; ++j)
          if (not finite(M(i, j))) throw std::invalid_argument(prefix + " of block " + name + " is not finite");
    }
  }

  g_w_t constant_gf(g_w_t const &F, block_matrices_t const &matrices) {
    g_w_t C = F;
    for (long b = 0; b < long(F.size()); ++b) {
      auto &d = C[b].data();
      for (long k = 0; k < d.extent(0); ++k)
        for (long i = 0; i < d.extent(1); ++i)
          for (long j = 0; j < d.extent(2); ++j) d(k, i, j) = matrices[b](i, j);
    }
    return C;
  }

} // namespace nrgljubljana_interface::dmft::detail
