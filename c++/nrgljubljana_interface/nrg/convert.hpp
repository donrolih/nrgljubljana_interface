// Conversions between TRIQS objects and NRG Ljubljana numerics types.
//
// NRG Ljubljana numerics works on one real function at a time (Table) and on spectral bins (SpectrumBins).
// The functions here take one part (real or imaginary) of one element of a refreq_log Green's function to a
// Table and back, and select one part of stored bins for broadening.

#pragma once

#include "./spectrum.hpp"

#include <nrgljubljana/numerics.hpp>

#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>

#include <vector>

namespace nrgljubljana_interface::nrg {

  /// Which part of a complex value a Table holds
  enum class part_t { real, imag };

  /// The mesh points, in ascending order (the Table x of a function on this mesh)
  std::vector<double> mesh_points(triqs::mesh::refreq_log const &mesh);

  /// Check that x are exactly the points of mesh. Throws std::invalid_argument naming the first mismatch.
  void check_on_mesh(std::vector<double> const &x, triqs::mesh::refreq_log const &mesh);

  /// One part of element (i, j) of g as a Table on the mesh points
  NRG::Numerics::Table to_table(triqs::gfs::gf_const_view<triqs::mesh::refreq_log, triqs::gfs::matrix_valued> g, long i, long j,
                                part_t part);

  /// One part of g as a Table on the mesh points
  NRG::Numerics::Table to_table(triqs::gfs::gf_const_view<triqs::mesh::refreq_log, triqs::gfs::scalar_valued> g, part_t part);

  /// Write t into one part of element (i, j) of g; the other part is kept. t must lie exactly on the mesh of g.
  void from_table(NRG::Numerics::Table const &t, triqs::gfs::gf_view<triqs::mesh::refreq_log, triqs::gfs::matrix_valued> g, long i, long j,
                  part_t part);

  /// Write t into one part of g; the other part is kept. t must lie exactly on the mesh of g.
  void from_table(NRG::Numerics::Table const &t, triqs::gfs::gf_view<triqs::mesh::refreq_log, triqs::gfs::scalar_valued> g, part_t part);

  /// The bins of one part of s, in the NRG Ljubljana type (input to broaden). part_t::imag requires a complex spectrum.
  NRG::Numerics::SpectrumBins to_nrg(spectrum_bins_t const &s, part_t part = part_t::real);

} // namespace nrgljubljana_interface::nrg
