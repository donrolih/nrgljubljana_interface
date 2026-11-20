// A lattice given by functions, e.g. written in Python.

#pragma once

#include "./lattice.hpp"

#include <functional>
#include <optional>
#include <string>
#include <utility>

namespace nrgljubljana_interface::dmft {

  /**
   * A lattice given by functions: the moments (G0, G1) of Σ and μ, optionally G0 alone, and a fixed on-site energy
   * ⟨h⟩ per block. The results must have the block structure and mesh of Σ and be finite; ⟨h⟩ must match the
   * structure it is asked for. The description is what resume compares, so it should identify the lattice. Errors
   * are std::invalid_argument.
   */
  class callable_lattice_t : public lattice_t {
   public:
    // Σ is passed by value, so a function written in Python gets its own copy
    using moments_t  = std::function<std::pair<g_w_t, g_w_t>(g_w_t Sigma, double mu)>;
    using local_gf_t = std::function<g_w_t(g_w_t Sigma, double mu)>;

    explicit callable_lattice_t(moments_t moments, block_matrices_t onsite_energy, std::optional<local_gf_t> local_gf = {},
                                std::string description = "callable");

    std::pair<g_w_t, g_w_t> moments(g_w_t const &Sigma, double mu) const override;
    g_w_t local_gf(g_w_t const &Sigma, double mu) const override;
    block_matrices_t onsite_energy(triqs::gfs::gf_struct_t const &structure) const override;
    [[nodiscard]] std::string description() const override { return description_; }

   private:
    moments_t moments_;
    block_matrices_t onsite_;
    std::optional<local_gf_t> local_gf_;
    std::string description_;
  };

} // namespace nrgljubljana_interface::dmft
