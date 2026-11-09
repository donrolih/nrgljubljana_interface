// Adapters between TRIQS objects and NRG Ljubljana numerics types (nrg/). No NRG run.

#include <gtest/gtest.h>

#include <nrgljubljana_interface/nrg/convert.hpp>
#include <nrgljubljana_interface/nrg/spectrum.hpp>

#include <nrgljubljana/numerics.hpp>

#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>

#include <complex>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace nrgljubljana_interface::nrg;
using triqs::gfs::gf;
using triqs::gfs::matrix_valued;
using triqs::gfs::scalar_valued;
using triqs::mesh::refreq_log;
using dcomplex = std::complex<double>;

namespace {

  refreq_log test_mesh() { return {1e-3, 1.0, 1.2}; }

  // Distinct values in each element and part, so a write to the wrong place shows
  dcomplex value(long k, long i, long j) { return {1.0 + k + 0.1 * i + 0.01 * j, -2.0 - k - 0.3 * i - 0.07 * j}; }

  gf<refreq_log, matrix_valued> matrix_gf() {
    auto g = gf<refreq_log, matrix_valued>{test_mesh(), {2, 2}};
    for (long k = 0; k < g.mesh().size(); ++k)
      for (long i = 0; i < 2; ++i)
        for (long j = 0; j < 2; ++j) g.data()(k, i, j) = value(k, i, j);
    return g;
  }

  gf<refreq_log, scalar_valued> scalar_gf() {
    auto g = gf<refreq_log, scalar_valued>{test_mesh()};
    for (long k = 0; k < g.mesh().size(); ++k) g.data()(k) = value(k, 0, 0);
    return g;
  }

  // A binary spectrum file in the working directory, removed at the end of the test
  struct binary_file {
    std::string path;
    binary_file(std::string name, std::vector<double> const &values) : path(std::move(name)) {
      std::ofstream out(path, std::ios::binary | std::ios::trunc);
      out.write(reinterpret_cast<char const *>(values.data()), static_cast<std::streamsize>(values.size() * sizeof(double)));
    }
    ~binary_file() { std::filesystem::remove(path); }
    binary_file(binary_file const &)            = delete;
    binary_file &operator=(binary_file const &) = delete;
  };

} // namespace

// ---- Mesh ----

TEST(NrgMesh, Points) {
  auto m = test_mesh();
  EXPECT_EQ(mesh_points(m), m.points());
  EXPECT_NO_THROW(check_on_mesh(mesh_points(m), m));
}

TEST(NrgMesh, RejectsOtherPoints) {
  auto m = test_mesh();
  auto x = mesh_points(m);
  x[3] += 1e-15;
  try {
    check_on_mesh(x, m);
    FAIL() << "A shifted point was accepted";
  } catch (std::invalid_argument const &e) { EXPECT_NE(std::string(e.what()).find("point 3"), std::string::npos) << e.what(); }

  x = mesh_points(m);
  x.pop_back();
  EXPECT_THROW(check_on_mesh(x, m), std::invalid_argument);
  EXPECT_THROW(check_on_mesh(mesh_points(refreq_log{1e-3, 1.0, 1.3}), m), std::invalid_argument);
}

// ---- Green's function element <-> Table ----

TEST(NrgTable, MatrixToTable) {
  auto g = matrix_gf();
  auto re = to_table(g, 1, 0, part_t::real);
  auto im = to_table(g, 1, 0, part_t::imag);
  EXPECT_EQ(re.x, g.mesh().points());
  EXPECT_EQ(im.x, g.mesh().points());
  for (long k = 0; k < g.mesh().size(); ++k) {
    EXPECT_EQ(re.y[k], value(k, 1, 0).real());
    EXPECT_EQ(im.y[k], value(k, 1, 0).imag());
  }
}

TEST(NrgTable, MatrixRoundTrip) {
  auto const g0 = matrix_gf();
  for (auto part : {part_t::real, part_t::imag}) {
    for (long i = 0; i < 2; ++i) {
      for (long j = 0; j < 2; ++j) {
        // Write a modified part into (i, j): only that part of that element changes
        auto t = to_table(g0, i, j, part);
        for (auto &y : t.y) y = 3 * y + 1;
        auto g = g0;
        from_table(t, g, i, j, part);
        for (long k = 0; k < g.mesh().size(); ++k) {
          for (long a = 0; a < 2; ++a) {
            for (long b = 0; b < 2; ++b) {
              auto expected = g0.data()(k, a, b);
              if (a == i and b == j) {
                if (part == part_t::real)
                  expected.real(3 * expected.real() + 1);
                else
                  expected.imag(3 * expected.imag() + 1);
              }
              EXPECT_EQ(g.data()(k, a, b), expected);
            }
          }
        }
        // And back: the Table read out equals the one written
        EXPECT_EQ(to_table(g, i, j, part).y, t.y);
      }
    }
  }
}

TEST(NrgTable, ScalarRoundTrip) {
  auto const g0 = scalar_gf();
  for (auto part : {part_t::real, part_t::imag}) {
    auto t = to_table(g0, part);
    EXPECT_EQ(t.x, g0.mesh().points());
    for (auto &y : t.y) y = -y;
    auto g = g0;
    from_table(t, g, part);
    for (long k = 0; k < g.mesh().size(); ++k) {
      auto expected = g0.data()(k);
      if (part == part_t::real)
        expected.real(-expected.real());
      else
        expected.imag(-expected.imag());
      EXPECT_EQ(g.data()(k), expected);
    }
    EXPECT_EQ(to_table(g, part).y, t.y);
  }
}

TEST(NrgTable, RejectsMismatch) {
  auto g = matrix_gf();
  auto const g0 = g;
  auto t = to_table(g, 0, 0, part_t::real);

  auto shifted = t;
  shifted.x[5] *= 1.0 + 1e-12;
  EXPECT_THROW(from_table(shifted, g, 0, 0, part_t::real), std::invalid_argument);

  auto short_y = t;
  short_y.y.pop_back();
  EXPECT_THROW(from_table(short_y, g, 0, 0, part_t::real), std::invalid_argument);

  EXPECT_THROW(from_table(t, g, 2, 0, part_t::real), std::out_of_range);
  EXPECT_THROW(to_table(g, 0, -1, part_t::real), std::out_of_range);

  auto s = scalar_gf();
  auto const s0 = s;
  EXPECT_THROW(from_table(shifted, s, part_t::imag), std::invalid_argument);

  // A rejected write leaves the Green's function unchanged
  EXPECT_EQ(g.data(), g0.data());
  EXPECT_EQ(s.data(), s0.data());
}

// ---- Spectrum bins ----

TEST(NrgSpectrum, LoadReal) {
  binary_file f("nrg_adapters_real.bin", {0.4, -2, -0.3, 1, 0.4, 0.25});
  auto s = load_spectrum(f.path, false);
  EXPECT_FALSE(s.is_complex);
  EXPECT_EQ(s.frequency, (std::vector<double>{0.4, -0.3, 0.4})); // file order, duplicates kept
  EXPECT_EQ(s.weight, (std::vector<double>{-2, 1, 0.25}));
  EXPECT_TRUE(s.weight_imag.empty());
}

TEST(NrgSpectrum, LoadComplex) {
  binary_file f("nrg_adapters_complex.bin", {0.4, -2, 3, -0.3, 1, -4, 0.4, 0.25, -0.5});
  auto s = load_spectrum(f.path, true);
  EXPECT_TRUE(s.is_complex);
  EXPECT_EQ(s.frequency, (std::vector<double>{0.4, -0.3, 0.4}));
  EXPECT_EQ(s.weight, (std::vector<double>{-2, 1, 0.25}));
  EXPECT_EQ(s.weight_imag, (std::vector<double>{3, -4, -0.5}));
}

TEST(NrgSpectrum, Empty) {
  binary_file f("nrg_adapters_empty.bin", {});
  for (bool is_complex : {false, true}) {
    auto s = load_spectrum(f.path, is_complex);
    EXPECT_EQ(s.is_complex, is_complex);
    EXPECT_TRUE(s.frequency.empty());
    EXPECT_TRUE(s.weight.empty());
    EXPECT_TRUE(s.weight_imag.empty());
    EXPECT_TRUE(to_nrg(s).empty());
  }
}

TEST(NrgSpectrum, ToNrgEqualsNrgLoader) {
  using NRG::Numerics::SpectrumLayout;
  auto same = [](NRG::Numerics::SpectrumBins const &a, NRG::Numerics::SpectrumBins const &b) {
    ASSERT_EQ(a.size(), b.size());
    for (std::size_t k = 0; k < a.size(); ++k) {
      EXPECT_EQ(a[k].frequency, b[k].frequency);
      EXPECT_EQ(a[k].weight, b[k].weight);
    }
  };

  binary_file r("nrg_adapters_real.bin", {0.4, -2, -0.3, 1, 0.4, 0.25});
  auto sr = load_spectrum(r.path, false);
  same(to_nrg(sr, part_t::real), NRG::Numerics::load_spectrum_bins(r.path, SpectrumLayout::real_pairs));
  EXPECT_THROW(to_nrg(sr, part_t::imag), std::invalid_argument);

  binary_file c("nrg_adapters_complex.bin", {0.4, -2, 3, -0.3, 1, -4, 0.4, 0.25, -0.5});
  auto sc = load_spectrum(c.path, true);
  same(to_nrg(sc, part_t::real), NRG::Numerics::load_spectrum_bins(c.path, SpectrumLayout::complex_real));
  same(to_nrg(sc, part_t::imag), NRG::Numerics::load_spectrum_bins(c.path, SpectrumLayout::complex_imaginary));

  // Inconsistent struct
  sc.weight_imag.pop_back();
  EXPECT_THROW(to_nrg(sc, part_t::imag), std::invalid_argument);
}

TEST(NrgSpectrum, Rejects) {
  EXPECT_THROW(load_spectrum("nrg_adapters_missing.bin", false), std::runtime_error);

  // Not a whole number of records
  binary_file odd("nrg_adapters_odd.bin", {0.1, 1, -1});
  EXPECT_THROW(load_spectrum(odd.path, false), std::invalid_argument);
  binary_file pair("nrg_adapters_pair.bin", {0.1, 1});
  EXPECT_THROW(load_spectrum(pair.path, true), std::invalid_argument);
}

TEST(NrgSpectrum, WrongLayoutNotAlwaysDetected) {
  // Two complex records are six doubles, which also read as three real records. The layout must come from the
  // caller; this documents that a wrong flag can pass unnoticed.
  binary_file f("nrg_adapters_complex2.bin", {0.4, -2, 3, -0.3, 1, -4});
  auto s = load_spectrum(f.path, false);
  EXPECT_EQ(s.frequency.size(), 3u);
}
