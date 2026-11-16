/*******************************************************************************
 *
 * nrgljubljana_interface: A TRIQS interface to the nrgljubliana impurity solver
 *
 * Copyright (c) 2019 The Simons foundation
 *   authors: Rok Zitko, Nils Wentzell
 *
 * nrgljubljana_interface is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * nrgljubljana_interface is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * nrgljubljana_interface. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#include "./solver_core.hpp"

#include <triqs/utility/exceptions.hpp> // TRIQS_RUNTIME_ERROR
#include <itertools/itertools.hpp>

#include <algorithm>  // max
#include <cmath>      // pow, log
#include <cstdlib>    // system
#include <sys/stat.h> // mkdir
#include <iostream>
#include <fstream>
#include <utility>
#include <deque>
#include <map>
#include <vector>
#include <iomanip>    // setprecision
#include <limits>     // max_digits10, MAX_PATH

#include <boost/lexical_cast.hpp>
#include <nrg-lib.hpp>

#include <mpi/mpi.hpp>
#include <mpi/string.hpp>

#include "openmp.hpp"
#include "templates/manifest.hpp"
#include "solver/param_file.hpp"

// std::filesystem is not production-ready as of 2020. As a workaround,
// we remove temporary files using "rm".
//
//#ifdef __has_include
//#  if __has_include(<filesystem>)
//#    include <filesystem>
//namespace fs = std::filesystem;
//#  elif __has_include(<experimental/filesystem>)
//#    include <experimental/filesystem>
//namespace fs = std::experimental::filesystem;
//#  else
//#    error "Missing <filesystem>"
//#  endif
//#endif

namespace nrgljubljana_interface {

  // Returns true if the GF structure is non-trivial, i.e., it has any elements at all.
  bool has_struct(const gf_struct_t &s) {
    unsigned int len = 0;
    for (const auto &[bl_name, bl_size] : s) len += bl_size;
    return len > 0;
  }

  namespace {

    std::string join(std::vector<std::string> const &items) {
      std::string result;
      for (auto const &item : items) result += (result.empty() ? "" : " ") + item;
      return result;
    }

  } // namespace

  void solver_core::read_manifest() {
    auto const dir = constr_params.get_model_dir();
    std::ifstream probe(dir + "/manifest");
    if (!probe) TRIQS_RUNTIME_ERROR << "Template " << dir << " has no manifest";

    templates::template_check_t check;
    try {
      manifest = templates::load_manifest(dir);
      check    = templates::check_against_template(manifest, dir);
    } catch (std::invalid_argument const &e) { TRIQS_RUNTIME_ERROR << e.what(); }
    correlator_type = check.correlator_type;
    auto const &m   = manifest;

    gf_struct    = m.gf;
    Delta_struct = m.delta;

    // chi_struct: the common block structure of the correlators other than G, F_l, F_r and I
    chi_struct.clear();
    std::string chi_name;
    for (auto const &[name, entries] : m.correlators) {
      if (templates::has_gf_structure(name)) continue;
      auto const s = m.correlator_structure(name);
      if (chi_name.empty()) {
        chi_struct = s;
        chi_name   = name;
      } else if (s != chi_struct)
        TRIQS_RUNTIME_ERROR << dir << "/manifest: [correlators] " << name << " and " << chi_name
                            << " have different block structures; chi_struct needs one";
    }

    // Spectrum lists of the NRG run, from the type of each correlator
    std::map<char, std::vector<std::string>> lists;
    for (auto const &[name, entries] : m.correlators)
      for (auto const &entry : entries) {
        auto &list = lists[check.correlator_type.at(entry.nrg_name)];
        if (std::find(list.begin(), list.end(), entry.nrg_name) == list.end()) list.push_back(entry.nrg_name);
      }
    constr_params.ops    = join(m.ops);
    constr_params.specs  = join(lists['s']);
    constr_params.specd  = join(lists['d']);
    constr_params.spect  = join(lists['t']);
    constr_params.specq  = join(lists['q']);
    constr_params.specot = join(lists['o']);
    constr_params.params = join(m.params);
    // Two coefficient channels per channel: a spin-polarized Wilson chain
    constr_params.polarized = m.layout == "uncoupled" and static_cast<long>(m.coefficient_channels.size()) == 2 * m.channels;
  }

  solver_core::solver_core(constr_params_t cp) : constr_params(cp) {
    if (world.rank() == 0)
      report_openMP(); // important for performance, so we report the settings here for easy inspection by the user

    read_manifest();

    // Create the hybridization function on a logarithmic mesh
    log_mesh = refreq_log{cp.mesh_min, cp.mesh_max, cp.mesh_ratio};
    // Note: Delta may have a different structure from G.
    Delta_w  = g_w_t{log_mesh, Delta_struct};
    // We also construct G_w and Sigma_w here to enable their initialization in the DMFT loops
    // prior to solver calls.
    if (has_struct(gf_struct)) {
      G_w = g_w_t{log_mesh, gf_struct};
      Sigma_w = g_w_t{log_mesh, gf_struct};
    }
  }

  // -------------------------------------------------------------------------------

  gf_struct_t solver_core::read_structure(const std::string &filename, bool mandatory = true) {
    std::ifstream F(constr_params.get_model_dir() + "/" + filename);
    if (mandatory && not F.is_open()) TRIQS_RUNTIME_ERROR << "Failed to open structure file " << filename;
    gf_struct_t _gf_struct;
    if (F) {
      std::string bl_name;
      int bl_size;
      while (F >> bl_name >> bl_size) { _gf_struct.emplace_back(bl_name, bl_size); }
    }
    return _gf_struct;
  }

  void solver_core::readGF(const std::string &name, std::optional<g_w_t> &G_w, gf_struct_t &_gf_struct) {
    if (!has_struct(_gf_struct))
      return;
    G_w = g_w_t{log_mesh, _gf_struct};
    for (int bl_idx : range(_gf_struct.size())) {
      long bl_size = Delta_w[bl_idx].target_shape()[0];
      for (auto [i, j] : product_range(bl_size, bl_size)) {
        auto bl_name     = Delta_w.block_names()[bl_idx];
        auto file_ending = std::string{"_"} + bl_name + "_" + std::to_string(i) + std::to_string(j) + ".dat";
        std::string imGfilename = "im" + name + file_ending;
        std::ifstream imG(imGfilename);
        std::string reGfilename = "re" + name + file_ending;
        std::ifstream reG(reGfilename);
        if (imG && reG) {
          double w, re, im;
          for (auto mp : log_mesh) {
            imG >> w >> im;
            reG >> w >> re;
            TRIQS_ASSERT2(abs(w-double(mp)) < 1e-8*abs(w), "frequency mismatch");
            (*G_w)[bl_idx][mp](i, j) = re + 1i * im;
          }
        } else {
          // If the files cannot be read (do not exist), the corresponding GF matrix component
          // is zero-ed out.
          for (auto mp : log_mesh) (*G_w)[bl_idx][mp](i, j) = 0.0;
        }
      }
    }
  }
  
  void solver_core::readA(const std::string &name, std::optional<g_w_t> &A_w, gf_struct_t &_gf_struct) {
    if (!has_struct(_gf_struct))
      return;
    A_w = g_w_t{log_mesh, _gf_struct};
    for (int bl_idx : range(_gf_struct.size())) {
      long bl_size = Delta_w[bl_idx].target_shape()[0];
      for (auto [i, j] : product_range(bl_size, bl_size)) {
        auto bl_name     = Delta_w.block_names()[bl_idx];
        auto file_ending = std::string{"_"} + bl_name + "_" + std::to_string(i) + std::to_string(j) + ".dat";
        std::string Afilename = name + file_ending;
        std::ifstream A(Afilename);
        if (A) {
          double w, re;
          for (auto mp : log_mesh) {
            A >> w >> re;
            TRIQS_ASSERT2(abs(w-double(mp)) < 1e-8*abs(w), "frequency mismatch");
            (*A_w)[bl_idx][mp](i, j) = re;
          }
        } else {
          for (auto mp : log_mesh) (*A_w)[bl_idx][mp](i, j) = 0.0;
        }
      }
    }
  }

  // Read expectation values and average over Nz runs
  void solver_core::readexpv(int Nz) {
    for (int cnt = 1; cnt <= Nz; cnt++) {
      const std::string expvfilename = std::to_string(cnt) + "/customfdm";
      std::ifstream F(expvfilename);
      if (!F) TRIQS_RUNTIME_ERROR << "Expectation values output file not found.";
      std::string snumber, skeyword, svalue;
      getline(F, snumber); // snumber not used
      getline(F, skeyword);
      getline(F, svalue);
      F.close();
      skeyword.erase(0,1); // drop #
      std::istringstream sskeyword(skeyword);
      std::deque<std::string> keywords{std::istream_iterator<std::string>{sskeyword},
        std::istream_iterator<std::string>{}};
      std::istringstream ssvalue(svalue);
      std::deque<double> values{std::istream_iterator<double>{ssvalue},
        std::istream_iterator<double>{}};
      TRIQS_ASSERT2(keywords.size() == values.size(), "customfdm corrupted");
      keywords.pop_front(); // ignore first column (temperature)
      values.pop_front();
      for (auto [k, v] : zip(keywords, values)) expv[k] += v; // zeroed by default constructor
    }
    for (auto &i : expv) { i.second /= Nz; } // calculate the average over discretization meshes
  }

  // Read thermodynamic variables (FDM algorithm) and average over Nz runs
  void solver_core::readtdfdm(int Nz) {
    for (int cnt = 1; cnt <= Nz; cnt++) {
      const std::string tdfdmfilename = std::to_string(cnt) + "/tdfdm";
      std::ifstream F(tdfdmfilename);
      if (!F) TRIQS_RUNTIME_ERROR << "Thermodynamic variables output file not found.";
      std::string skeyword, svalue;
      getline(F, skeyword);
      getline(F, svalue);
      F.close();
      skeyword.erase(0,1); // drop #
      std::istringstream sskeyword(skeyword);
      std::deque<std::string> keywords{std::istream_iterator<std::string>{sskeyword},
        std::istream_iterator<std::string>{}};
      std::istringstream ssvalue(svalue);
      std::deque<double> values{std::istream_iterator<double>{ssvalue},
        std::istream_iterator<double>{}};
      TRIQS_ASSERT2(keywords.size() == values.size(), "tdfdm corrupted");
      keywords.pop_front(); // ignore first column (temperature)
      values.pop_front();
      for (auto [k, v] : zip(keywords, values)) tdfdm[k] += v; // zeroed by default constructor
    }
    for (auto &i : tdfdm) { i.second /= Nz; } // calculate the average over discretization meshes
  }

  inline void call(const std::string &command, bool verbose = true, bool exit_on_failure = true) {
    std::string s = command + (verbose ? "" : " >/dev/null");
    std::cout << "running " << s << std::endl;
    if (system(s.c_str()) != 0) {
      if (exit_on_failure) 
        TRIQS_RUNTIME_ERROR << "Failure running NRGLjubljana_interface script: " << s;
      else
        std::cerr << "Warning: failure running NRGLjubljana_interface script: " << s << std::endl;
    }
  }

  // Write Gamma=-Im(Delta_w) to a file
  void solver_core::write_gamma() {
    for (int bl_idx : range(Delta_struct.size())) {
      long bl_size = Delta_w[bl_idx].target_shape()[0];
      auto bl_name = Delta_w.block_names()[bl_idx];
      for (auto [i, j] : product_range(bl_size, bl_size)) {
        std::ofstream F("Gamma_" + bl_name + "_" + std::to_string(i) + std::to_string(j) + ".dat");
        F << std::setprecision(std::numeric_limits<double>::max_digits10); // number of decimal digits necessary to differentiate all values of this type
        const double cutoff = 1e-8; // ensure hybridisation function is positive
        for (auto w : Delta_w[bl_idx].mesh()) {
          double value = -Delta_w[bl_idx][w](i, j).imag();
          if (value < cutoff) { value = cutoff; }
          F << double(w) << " " << value << std::endl;
        }
      }
    }
  }

  inline std::string get_current_dir() {
    char temp[PATH_MAX]; // on stack, no memory leak
    return ( getcwd(temp, sizeof(temp)) ? std::string(temp) : std::string("") );
  }
   
  void solver_core::solve(solve_params_t const &sp) {
    request = solve_request_t{};
    try {
      validate_request(request, manifest);
      if (check_method(sp, request) and world.rank() == 0)
        std::cerr << "Warning: solve_params.method is deprecated; FDM is the only method and always used\n";
    } catch (std::invalid_argument const &e) { TRIQS_RUNTIME_ERROR << e.what(); }
    auto const z = twist_values(request, sp.Nz);
    last_solve_params = sp;
    container_set::operator=(container_set{});  // Reset the results
    std::string tempdir{};
    if (world.rank() == 0) {
      std::cout << "Starting NRG solver\n";
      tempdir = create_tempdir("");  // Create a temporary directory for file-based interfacing (RAM disk is best)
    }
    world.barrier(); // Ensures temporary directory is created
    mpi::broadcast(tempdir, world);
    std::string cwd = get_current_dir();
    if (chdir(tempdir.c_str()) != 0) TRIQS_RUNTIME_ERROR << "chdir to tempdir failed.";
    if (world.rank() == 0) {
      write_gamma();
      // we need a mock param file for 'adapt' tool
      generate_param_file(1.0);
      call(constr_params.get_model_dir() + "/prepare", verbose);
      call("./discretize", verbose);
      for (auto k : range(z.size())) instantiate(z[k], std::to_string(k + 1));
    }
    world.barrier(); // Ensures the initialization scripts have completed

    // Perform the calculations in parallel
    auto chunk = mpi::chunk(range(1, long(z.size()) + 1), world);
    for (int i: chunk)
      solve_one(std::to_string(i));

    world.barrier(); // Ensures all subcalculations are completed
    // Post-Processing in perl script
    if (world.rank() == 0) {
      std::cout << "Post-processing results" << std::endl;
      call("./process", verbose);
    }

    world.barrier(); // Ensures post-processing is completed
    // Read Results into TRIQS Containers
    if (world.rank() == 0)
      std::cout << "Reading results from files" << std::endl;
    readexpv(int(z.size()));
    readtdfdm(int(z.size()));
    readGF("G", G_w, gf_struct);
    readGF("F_l", F_l_w, gf_struct);
    readGF("F_r", F_r_w, gf_struct);
    readGF("I", I_w, gf_struct);
    readGF("SigmaHartree", SigmaHartree_w, gf_struct);
    readA("A", A_w, gf_struct);
    readA("B_l", B_l_w, gf_struct);
    readA("B_r", B_r_w, gf_struct);
    readA("C", C_w, gf_struct);
    readGF("SS", chi_SS_w, chi_struct);
    readGF("NN", chi_NN_w, chi_struct);

    // Post-Processing in C++ interface
    if (has_struct(gf_struct)) {
      Sigma_w = (*SigmaHartree_w) + (*I_w) - (*F_l_w) / (*G_w) * (*F_r_w);
    }
    // Cleanup
    world.barrier(); // Ensures all processes have read the results before cleanup
    if (chdir(cwd.c_str()) != 0) TRIQS_RUNTIME_ERROR << "failed to return from tempdir";
    world.barrier(); // Ensures all processes have chdired from temp dir before it is deleted
    if (world.rank() == 0 && !keep_temp_dir) {
        //      std::error_code ec;
        //      fs::remove_all(tempdir, ec);
        //      if (ec) std::cerr << "Warning: failed to remove the temporary directory." << std::endl;
        //    Workaround:
        call("rm -rf " + tempdir, verbose, false);
    }
  }

  // -------------------------------------------------------------------------------

  void solver_core::set_params() {
    const constr_params_t &cp = constr_params;
    solve_params_t &sp        = *last_solve_params;
    nrg_params_t &np          = nrg_params; // only nrg_params allowed to be changed!

    // Test if the low-level paramerers are sensible for use with the
    // high-level interface.
    if (np.discretization != "Z") TRIQS_RUNTIME_ERROR << "Must use discretization=Z in the high-level solver interface.";

    // Automatically set (override) some low-level parameters
    if (sp.Tmin > 0) { // If Tmin is set, determine the required length of the Wilson chain.
      np.Nmax    = 0;
      auto scale = [=](int n) {
        return (1. - 1. / sp.Lambda) / std::log(sp.Lambda) * std::pow(sp.Lambda, -(np.z - 1)) * std::pow(sp.Lambda, -(n - 1) / 2.);
      };
      while (scale(np.Nmax + 1) >= sp.Tmin) np.Nmax++;
    }
    if (np.mMAX < 0) // Determine the number of sites in the star representation
      np.mMAX = std::max(80, 2 * np.Nmax);
    if (np.xmax < 0) // Length of the x-interval in the discretization=Z (ODE) approach
      np.xmax = np.Nmax / 2. + 2.;
    if (sp.bandrescale < 0) // Make the NRG energy window correspond to the extent of the frequency mesh
      sp.bandrescale = cp.mesh_max;
  }

  void solver_core::set_nrg_params(nrg_params_t const &nrg_params_) { nrg_params = nrg_params_; }

  // Creates input files for NRG iteration (param & data) by calling the "instantiate" script
  void solver_core::instantiate(double z, const std::string &taskdir) {
    std::cout << "Preparing z=" << z << " taskdir=" << taskdir << std::endl;
    TRIQS_ASSERT(world.rank() == 0); // should not run in parallel!
    generate_param_file(z);
    if (mkdir(taskdir.c_str(), 0755) != 0) TRIQS_RUNTIME_ERROR << "failed to mkdir taskdir " << taskdir;
    call("./instantiate " + taskdir, verbose);
  }

  // Solve the problem for a single value of the twist parameter z
  void solver_core::solve_one(const std::string &taskdir) {
    std::cout << "Solving taskdir=" << taskdir << " rank=" << world.rank() << std::endl;
    if (chdir(taskdir.c_str()) != 0) TRIQS_RUNTIME_ERROR << "failed to chdir to taskdir " << taskdir;
    // Solve the impurity model
    std::streambuf *old = std::cout.rdbuf();
    std::ofstream log_stream;
    if (!verbose) {
      log_stream.open("log");
      if (!log_stream) TRIQS_RUNTIME_ERROR << "failed to open log file";
      std::cout.rdbuf(log_stream.rdbuf()); // redirect stdout to log file
    }
    NRG::set_workdir(""); // Defaults to . but may be overridden by NRG_WORKDIR in the environment
    NRG::run_nrg_master("");
    if (!verbose) std::cout.rdbuf(old); // restore
    if (chdir("..") != 0) TRIQS_RUNTIME_ERROR << "failed to return from taskdir " << taskdir;
  }

  std::string solver_core::create_tempdir(const std::string &tempdir_ = "") {
    const std::string default_tempdir = ".";
    std::string tempdir = default_tempdir;
    if (const char *env_w = std::getenv("NRG_TEMPDIR")) tempdir = env_w;
    if (!tempdir_.empty()) tempdir = tempdir_;
    const std::string tempdir_template = tempdir + "/nrg_tempdir_XXXXXX";
    size_t len = tempdir_template.length()+1;
    auto x = std::make_unique<char[]>(len); // NOLINT
    strncpy(x.get(), tempdir_template.c_str(), len);
    if (auto w = mkdtemp(x.get())) // create a unique directory
      return w;
    else
      TRIQS_RUNTIME_ERROR << "Failed to create a directory for temporary files.";
  }
   
  void solver_core::check_model_params(const solve_params_t &sp) {
    std::istringstream ss(constr_params.params);
    do {
      std::string p;
      ss >> p;
      if (p != "" && sp.model_parameters.count(p) != 1)
        TRIQS_RUNTIME_ERROR << "Model parameter " << p << " has not been defined.";
    } while (ss);
  }

  void solver_core::generate_param_file(double z) {
    const constr_params_t &cp = constr_params;
    const solve_params_t &sp  = *last_solve_params;
    nrg_params_t &np          = nrg_params;
    np.z                      = z;
    // Automatically establish appropriate default values for the high-level interface
    set_params();
    // Check if all model parameters have been defined
    check_model_params(sp);
    // Generate the parameter file; the request decides the spectrum lists and the density-matrix settings
    auto cp_solve = cp;
    apply_request(request, manifest, correlator_type, cp_solve, np);
    std::ofstream F("param");
    F << render_param_file(cp_solve, sp, np);
  }

  // -------------------------------------------------------------------------------

  // Function that writes a solver object to hdf5 file

  void h5_write(h5::group h5group, std::string subgroup_name, solver_core const &s) {
    auto grp = h5group.create_group(subgroup_name);
    h5_write_attribute(grp, "Format", solver_core::hdf5_format());
    h5_write_attribute(grp, "TRIQS_GIT_HASH", std::string(AS_STRING(TRIQS_GIT_HASH)));
    h5_write_attribute(grp, "NRGLJUBLJANA_INTERFACE_GIT_HASH", std::string(AS_STRING(NRGLJUBLJANA_INTERFACE_GIT_HASH)));
    h5_write(grp, "", s.result_set());
    h5_write(grp, "constr_params", s.constr_params);
    h5_write(grp, "nrg_params", s.nrg_params);
    h5_write(grp, "last_solve_params", s.last_solve_params);
    h5_write(grp, "Delta_w", s.Delta_w); // !!
  }

  // Function that constructs a solver object from an hdf5 file
  solver_core solver_core::h5_read_construct(h5::group h5group, std::string subgroup_name) {
    auto grp           = h5group.open_group(subgroup_name);
    auto constr_params = h5_read<constr_params_t>(grp, "constr_params");
    auto s             = solver_core{constr_params};
    h5_read(grp, "", s.result_set());
    h5_read(grp, "nrg_params", s.nrg_params);
    h5_read(grp, "last_solve_params", s.last_solve_params);
    h5_read(grp, "Delta_w", s.Delta_w);
    return s;
  }

  // Hilbert transform for refreq objects
  std::complex<double> hilbert_transform_refreq(const c_w_cvt &gf, std::complex<double> z){
    return hilbert_transform(gf, z);
  }

  matrix<std::complex<double>> hilbert_transform_elementwise(const m_w_cvt &gf, std::complex<double> z) {
    return triqs::gfs::hilbert_transform_elementwise(gf, z);
  }

} // namespace nrgljubljana_interface
