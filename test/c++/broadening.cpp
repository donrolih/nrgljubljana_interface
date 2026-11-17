// Broadening and Kramers-Kronig through NRG Ljubljana numerics (nrg/broadening.hpp) agree with the broaden and kk
// tools, also for spectra with negative weights (off-diagonal correlators). Kramers-Kronig is the analytic kernel on
// Steffen interpolation, which the kk tool computes with "-a analytic -i steffen" (its default is QAG on Akima).

#include <gtest/gtest.h>

#include <nrgljubljana_interface/nrg/broadening.hpp>
#include <nrgljubljana_interface/nrg/output.hpp>

#include <triqs/mesh.hpp>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace nrgljubljana_interface;
namespace fs = std::filesystem;

namespace {

  // Bins of two twists; the second has negative weights, as n_d1-n_d2 in 2orb-UJ
  std::vector<nrg::spectrum_bins_t> twist_bins() {
    nrg::spectrum_bins_t a, b;
    for (int k = -40; k <= 40; ++k) {
      if (k == 0) continue;
      double const w = std::copysign(1e-4 * std::pow(1.25, std::abs(k)), double(k));
      a.frequency.push_back(w);
      a.weight.push_back(0.01 / (1 + w * w));
      b.frequency.push_back(w * 1.07);
      b.weight.push_back(w > 0 ? -0.05 : 0.006 * std::exp(-std::abs(w)));
    }
    return {a, b};
  }

  void write_bins(std::string const &file, nrg::spectrum_bins_t const &s) {
    std::ofstream f(file, std::ios::binary);
    for (std::size_t k = 0; k < s.frequency.size(); ++k) {
      f.write(reinterpret_cast<char const *>(&s.frequency[k]), sizeof(double));
      f.write(reinterpret_cast<char const *>(&s.weight[k]), sizeof(double));
    }
  }

  // Two-column text file of a tool: x and y
  std::pair<std::vector<double>, std::vector<double>> read_columns(std::string const &file) {
    std::istringstream in(nrg::read_text_file(file));
    std::vector<double> x, y;
    for (double a = 0, b = 0; in >> a >> b;) {
      x.push_back(a);
      y.push_back(b);
    }
    return {x, y};
  }

  void run(std::string const &command) { ASSERT_EQ(std::system(command.c_str()), 0) << command; }

  struct workdir_t {
    fs::path old = fs::current_path();
    fs::path dir = fs::temp_directory_path() / ("broadening_test_" + std::to_string(::getpid()));
    workdir_t() {
      fs::remove_all(dir);
      fs::create_directories(dir);
      fs::current_path(dir);
    }
    ~workdir_t() {
      fs::current_path(old);
      fs::remove_all(dir);
    }
  };

  double const mesh_min = 1e-3, mesh_max = 1.0, mesh_ratio = 1.1;

} // namespace

TEST(Broadening, SameAsBroadenTool) {
  workdir_t const work;
  auto const twists = twist_bins();
  for (int i = 0; i < 2; ++i) {
    fs::create_directory(std::to_string(i + 1));
    write_bins(std::to_string(i + 1) + "/spec.bin", twists[i]);
  }
  double const alpha = 0.4, T = 1e-3, gamma = 0.2;
  std::ostringstream cmd;
  cmd << std::setprecision(17) << "broaden -x " << gamma << " -m " << mesh_min << " -M " << mesh_max << " -r " << mesh_ratio << " spec.bin 2 "
      << alpha << " " << T << " 1e-9 > broaden.log 2>&1";
  run(cmd.str());
  auto const [x, tool] = read_columns("spec.dat");

  auto const points = triqs::mesh::refreq_log{mesh_min, mesh_max, mesh_ratio}.points();
  ASSERT_EQ(x.size(), points.size());
  for (std::size_t k = 0; k < x.size(); ++k) ASSERT_NEAR(x[k], points[k], 1e-15 * std::abs(points[k])) << k;

  auto const ours = nrg::broaden({&twists[0], &twists[1]}, points, {alpha, T, gamma});
  bool negative = false;
  for (std::size_t k = 0; k < ours.size(); ++k) {
    EXPECT_NEAR(ours[k], tool[k], 1e-14 * std::max(1.0, std::abs(tool[k]))) << "w = " << points[k];
    negative = negative or ours[k] < 0;
  }
  EXPECT_TRUE(negative) << "the test spectrum should have negative parts";
}

TEST(Broadening, SameAsKkTool) {
  workdir_t const work;
  auto const points = triqs::mesh::refreq_log{mesh_min, mesh_max, mesh_ratio}.points();
  std::vector<double> im(points.size());
  for (std::size_t k = 0; k < points.size(); ++k) im[k] = -0.1 / (points[k] * points[k] + 0.01) + (points[k] > 0.3 ? 0.02 : 0.0);
  {
    std::ofstream f("im.dat");
    f << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (std::size_t k = 0; k < points.size(); ++k) f << points[k] << " " << im[k] << "\n";
  }
  run("kk -a analytic -i steffen im.dat re.dat > kk.log 2>&1"); // the method of NRG Ljubljana numerics
  auto const [x, tool] = read_columns("re.dat");
  ASSERT_EQ(x.size(), points.size());

  auto const ours = nrg::kramers_kronig(points, im);
  for (std::size_t k = 0; k < ours.size(); ++k) EXPECT_NEAR(ours[k], tool[k], 1e-14 * std::max(1.0, std::abs(tool[k]))) << "w = " << points[k];
}
