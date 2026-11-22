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
#include "./dmft/causality.hpp"
#include "./nrg/chain.hpp"
#include "./nrg/parameter_file.hpp"

#include <triqs/utility/exceptions.hpp> // TRIQS_RUNTIME_ERROR
#include <itertools/itertools.hpp>

#include <algorithm>  // max
#include <cmath>      // pow, log
#include <cstdlib>    // system
#include <filesystem>
#include <sstream>
#include <sys/stat.h> // mkdir
#include <iostream>
#include <fstream>
#include <utility>
#include <deque>
#include <map>
#include <numbers>
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
#include "solver/collect.hpp"
#include "nrg/revision.hpp"

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

  solver_core::solver_core(constr_params_t cp, bool mpi)
     : solver_core(std::move(cp), mpi ? ::mpi::communicator{} : ::mpi::communicator{MPI_COMM_SELF}) {}

  solver_core::solver_core(constr_params_t cp, mpi::communicator c) : comm(c), constr_params(std::move(cp)) {
    if (comm.rank() == 0)
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

  inline void call(const std::string &command, bool verbose = true) {
    std::string s = command + (verbose ? "" : " >/dev/null");
    std::cout << "running " << s << std::endl;
    if (system(s.c_str()) != 0) TRIQS_RUNTIME_ERROR << "Failure running NRGLjubljana_interface script: " << s;
  }

  // Write Gamma=-Im(Delta_w) to a file
  void solver_core::write_gamma() {
    for (int bl_idx : range(Delta_struct.size())) {
      long bl_size = Delta_w[bl_idx].target_shape()[0];
      auto bl_name = Delta_w.block_names()[bl_idx];
      for (auto [i, j] : product_range(bl_size, bl_size)) {
        std::ofstream F("Gamma_" + bl_name + "_" + std::to_string(i) + std::to_string(j) + ".dat");
        F << std::setprecision(std::numeric_limits<double>::max_digits10); // number of decimal digits necessary to differentiate all values of this type
        const double cutoff = last_solve_params->Gamma_floor; // ensure hybridisation function is positive
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
   
  void solver_core::fill_containers(result_t const &r) {
    expv  = r.expv;
    tdfdm = r.tdfdm;
    // From the result, or zero if the template does not provide it
    auto const set = [&](std::optional<g_w_t> &member, std::string const &name, gf_struct_t const &s, bool provided) {
      member.reset();
      if (not has_struct(s)) return;
      if (auto it = r.processed.find(name); it != r.processed.end()) {
        member = it->second;
      } else if (not provided) {
        member = g_w_t{log_mesh, s};
        for (auto &b : member->data()) b.data() = 0;
      }
    };
    auto const provides = [&](std::string const &name) { return manifest.correlators.contains(name); };
    bool const self_energy = std::ranges::all_of(templates::gf_correlators, provides);
    set(G_w, "G", gf_struct, provides("G"));
    set(F_l_w, "F_l", gf_struct, provides("F_l"));
    set(F_r_w, "F_r", gf_struct, provides("F_r"));
    set(I_w, "I", gf_struct, provides("I"));
    set(SigmaHartree_w, "SigmaHartree", gf_struct, self_energy);
    set(chi_NN_w, "chi_NN", chi_struct, provides("chi_NN"));
    set(chi_SS_w, "chi_SS", chi_struct, provides("chi_SS"));

    // Spectral densities -Im X / pi of G, F_l, F_r and I
    auto const density = [](std::optional<g_w_t> &member, std::optional<g_w_t> const &g) {
      member.reset();
      if (not g) return;
      member = *g;
      for (auto &b : member->data())
        for (auto &x : b.data()) x = -x.imag() / std::numbers::pi;
    };
    density(A_w, G_w);
    density(B_l_w, F_l_w);
    density(B_r_w, F_r_w);
    density(C_w, I_w);

    Sigma_w.reset();
    if (auto it = r.processed.find("Sigma"); it != r.processed.end()) Sigma_w = it->second;
  }

  g_w_t solver_core::broaden(std::string const &correlator, broadening_params_t const &b) {
    if (not result) TRIQS_RUNTIME_ERROR << "broaden: no result; solve first";
    try {
      result->processed[correlator] = nrgljubljana_interface::broaden(*result, correlator, b, log_mesh);
    } catch (std::exception const &e) { TRIQS_RUNTIME_ERROR << "broaden: " << e.what(); }
    fill_containers(*result);
    return result->processed.at(correlator);
  }

  void solver_core::sync() {
    if (comm.size() > 1) comm.barrier();
  }

  void solver_core::solve(solve_params_t const &sp) { solve(sp, options_t{}); }

  void solver_core::solve(solve_params_t const &sp, options_t const &opts) {
    // nrg_params as given before and after the solve: the values resolved for each twist are in its param file
    struct restore_t {
      nrg_params_t &target;
      nrg_params_t saved;
      ~restore_t() { target = saved; }
    } restore_nrg_params{nrg_params, nrg_params};

    options = opts;
    try {
      validate_options(options, manifest);
      if (check_method(sp, options) and comm.rank() == 0)
        std::cerr << "Warning: solve_params.method is deprecated; FDM is the only method and always used\n";
    } catch (std::invalid_argument const &e) { TRIQS_RUNTIME_ERROR << e.what(); }
    auto const z = twist_values(options, sp.Nz);
    last_solve_params = sp;
    if (last_solve_params->bandrescale < 0) last_solve_params->bandrescale = constr_params.mesh_max; // on every rank
    result.reset();

    // What was solved, as given (the values resolved for each twist are in its param file)
    result_t r;
    r.options            = options;
    r.solve_params       = sp;
    r.nrg_params         = nrg_params;
    r.model              = constr_params.model;
    r.symtype            = constr_params.symtype;
    r.manifest           = manifest.text;
    for (auto const &[corr, t] : correlator_type) r.correlator_type[corr] = std::string(1, t);
    r.nrg_revision       = nrg::revision();
    r.mesh               = log_mesh;
    r.z                  = z;
    r.complete_z_average = complete_z_average(z, sp.Nz);
    container_set::operator=(container_set{});  // Reset the results
    std::string tempdir{};
    if (comm.rank() == 0) {
      std::cout << "Starting NRG solver\n";
      tempdir = std::filesystem::absolute(create_tempdir("")).string(); // file-based interfacing (RAM disk is best)
    }
    sync(); // Ensures temporary directory is created
    if (comm.size() > 1) mpi::broadcast(tempdir, comm);
    std::string cwd    = get_current_dir();
    auto const spectra = selected_spectra(options, manifest, correlator_type);

    if (manifest.pipeline == "library") {
      solve_library(z, tempdir, spectra, r);
    } else {
      if (chdir(tempdir.c_str()) != 0) TRIQS_RUNTIME_ERROR << "chdir to tempdir failed.";
      if (comm.rank() == 0) {
        write_gamma();
        // we need a mock param file for 'adapt' tool
        generate_param_file(1.0);
        call(constr_params.get_model_dir() + "/prepare", verbose);
        call("./discretize", verbose);
        for (auto k : range(z.size())) instantiate(z[k], std::to_string(k + 1));
      }
      sync(); // Ensures the initialization scripts have completed

      // Perform the calculations in parallel
      auto chunk = mpi::chunk(range(1, long(z.size()) + 1), comm);
      for (int i : chunk) solve_one(std::to_string(i));

      sync(); // Ensures all subcalculations are completed

      // Raw results: every rank reads all twists (the working directory is shared)
      for (auto k : range(z.size())) r.twists.push_back(read_twist(std::to_string(k + 1), z[k], spectra, options.fdm, options.energies > 0));
    }

    // Spectra and averages, on every rank from the same input
    try {
      process(r, manifest);
    } catch (std::exception const &e) { TRIQS_RUNTIME_ERROR << "Processing the results failed: " << e.what(); }
    if (not options.keep_bins)
      for (auto &t : r.twists) t.spectra.clear();
    fill_containers(r);
    result = std::move(r);

    // Cleanup
    sync(); // Ensures all processes have read the results before cleanup
    if (chdir(cwd.c_str()) != 0) TRIQS_RUNTIME_ERROR << "failed to return from tempdir";
    sync(); // Ensures all processes have chdired from temp dir before it is deleted
    if (comm.rank() == 0 && !keep_temp_dir && !options.keep_workdir) {
      std::error_code ec;
      std::filesystem::remove_all(tempdir, ec);
      if (ec) std::cerr << "Warning: failed to remove the temporary directory " << tempdir << ": " << ec.message() << std::endl;
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
    TRIQS_ASSERT(comm.rank() == 0); // should not run in parallel!
    generate_param_file(z);
    if (mkdir(taskdir.c_str(), 0755) != 0) TRIQS_RUNTIME_ERROR << "failed to mkdir taskdir " << taskdir;
    call("./instantiate " + taskdir, verbose);
  }

  // Solve the problem for a single value of the twist parameter z
  void solver_core::solve_one(const std::string &taskdir) {
    std::cout << "Solving taskdir=" << taskdir << " rank=" << comm.rank() << std::endl;
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

  std::string solver_core::param_text(double z) {
    const constr_params_t &cp = constr_params;
    const solve_params_t &sp  = *last_solve_params;
    nrg_params_t &np          = nrg_params;
    np.z                      = z;
    // Automatically establish appropriate default values for the high-level interface
    set_params();
    // Check if all model parameters have been defined
    check_model_params(sp);
    // The options decide the spectrum lists and the density-matrix settings
    auto cp_solve = cp;
    apply_options(options, manifest, correlator_type, cp_solve, np);
    return render_param_file(cp_solve, sp, np);
  }

  void solver_core::generate_param_file(double z) {
    std::ofstream F("param");
    F << param_text(z);
  }

  namespace {

    // Run f on every rank; an exception on any rank is raised on every rank (the message of the lowest failing rank),
    // so that none is left waiting in a collective call
    template <typename F> void collectively(mpi::communicator const &comm, F &&f) {
      std::string error;
      try {
        f();
      } catch (std::exception const &e) { error = e.what(); }
      if (comm.size() > 1) {
        int const first = mpi::all_reduce(error.empty() ? comm.size() : comm.rank(), comm, MPI_MIN);
        if (first < comm.size()) mpi::broadcast(error, comm, first);
      }
      if (not error.empty()) TRIQS_RUNTIME_ERROR << error;
    }

    std::string number(double x) {
      std::ostringstream os;
      os.precision(17);
      os << x;
      return os.str();
    }

    // The text of parameter-file sections
    std::string section_text(nrg::parameter_sections_t const &sections) {
      std::string text;
      for (auto const &[name, keys] : sections) {
        text += "[" + name + "]\n";
        for (auto const &[key, value] : keys) text += key + "=" + value + "\n";
      }
      return text;
    }

    // The coefficients of all twists from rank 0 to every rank
    void broadcast(std::vector<nrg::chain_coefficients_t> &chains, mpi::communicator const &comm) {
      if (comm.size() == 1) return;
      long n = long(chains.size());
      mpi::broadcast(n, comm, 0);
      chains.resize(n);
      for (auto &c : chains) {
        mpi::broadcast(c.z, comm, 0);
        mpi::broadcast(c.Lambda, comm, 0);
        mpi::broadcast(c.bandrescale, comm, 0);
        long channels = long(c.channels.size());
        mpi::broadcast(channels, comm, 0);
        c.channels.resize(channels);
        for (auto &ch : c.channels) {
          mpi::broadcast(ch.theta, comm, 0);
          mpi::broadcast(ch.xi, comm, 0);
          mpi::broadcast(ch.zeta, comm, 0);
        }
      }
    }

  } // namespace

  void solver_core::solve_library(std::vector<double> const &z, std::string const &work, std::map<std::string, std::string> const &spectra,
                                  result_t &r) {
    namespace fs              = std::filesystem;
    solve_params_t const &sp  = *last_solve_params;
    std::string const pre     = "in-process solve: ";

    // Step 2: Γ = (Δ† − Δ)/2i, its eigenvalues raised to the floor
    if (not(std::isfinite(sp.Gamma_floor) and sp.Gamma_floor >= 0)) TRIQS_RUNTIME_ERROR << pre << "Gamma_floor must be finite and >= 0";
    auto const [Delta_p, projection] = dmft::causal_projection(Delta_w, sp.Gamma_floor);
    r.Gamma_projected_points          = projection.projected_points;
    r.Gamma_largest_violation         = projection.largest_violation;
    auto const Gamma                  = dmft::spectral_part(Delta_p);

    // Step 3: Γ of the coefficient channels; any other nonzero element is an error, never dropped
    std::vector<std::pair<long, long>> channel; // (block, orbital) of each coefficient channel
    for (auto const &e : manifest.coefficient_channels) {
      if (e.i != e.j) TRIQS_RUNTIME_ERROR << pre << "the coefficient channel " << templates::to_string(e) << " is off-diagonal";
      auto const &names = Gamma.block_names();
      auto const b      = std::find(names.begin(), names.end(), e.block) - names.begin();
      if (b == long(names.size())) TRIQS_RUNTIME_ERROR << pre << "no block " << e.block << " in Δ";
      channel.emplace_back(b, e.i);
    }
    for (long b = 0; b < long(Gamma.size()); ++b) {
      auto const &d = Gamma[b].data();
      for (long i = 0; i < d.extent(1); ++i)
        for (long j = 0; j < d.extent(2); ++j) {
          bool const is_channel = i == j and std::find(channel.begin(), channel.end(), std::pair<long, long>{b, i}) != channel.end();
          if (is_channel) continue;
          double largest = 0;
          for (long k = 0; k < d.extent(0); ++k) largest = std::max(largest, std::abs(d(k, i, j)));
          if (largest > 0)
            TRIQS_RUNTIME_ERROR << pre << "Γ of " << Gamma.block_names()[b] << "[" << i << "," << j << "] is nonzero (largest " << largest
                                << "), which the uncoupled channels of the template cannot represent";
        }
    }
    nrg::channel_gamma_t gamma;
    {
      auto const &pts = Gamma[0].mesh().points();
      gamma.omega.assign(pts.begin(), pts.end());
      gamma.gamma.assign(pts.size(), std::vector<double>(channel.size()));
      for (std::size_t c = 0; c < channel.size(); ++c) {
        auto const &[b, i] = channel[c];
        for (std::size_t k = 0; k < pts.size(); ++k) gamma.gamma[k][c] = Gamma[b].data()(long(k), i, i).real();
      }
    }

    // Step 4: the parameter sections of each twist; Nmax from Tmin (with bandrescale) unless Tmin <= 0
    std::vector<nrg::parameter_sections_t> sections(z.size());
    std::vector<unsigned> nmax(z.size());
    collectively(comm, [&] {
      for (std::size_t k = 0; k < z.size(); ++k) {
        sections[k] = nrg::read_parameter_sections(param_text(z[k]), "param");
        auto &param = sections[k].at("param");
        param["density_interpolation"] = nrg_params.density_interpolation;
        if (sp.Tmin > 0) {
          auto with_tmin = param;
          with_tmin.erase("Nmax");
          with_tmin["Tmin"] = number(sp.Tmin);
          nmax[k]           = nrg::chain_length(with_tmin);
        } else {
          if (nrg_params.Nmax < 1) throw std::invalid_argument("with Tmin <= 0 the chain length needs Nmax (set_nrg_params) >= 1");
          nmax[k] = unsigned(nrg_params.Nmax);
        }
        param["Nmax"] = std::to_string(nmax[k]);
      }
      unsigned const longest = *std::max_element(nmax.begin(), nmax.end());
      for (auto &s : sections) {
        s.at("param")["mMAX"] = std::to_string(std::max(80u, 2 * longest));
        nrg::validate_chain_parameters(s.at("param"));
      }
    });
    unsigned const longest = *std::max_element(nmax.begin(), nmax.end());

    // Step 5: the chains of all twists on rank 0
    std::vector<nrg::chain_coefficients_t> chains;
    collectively(comm, [&] {
      if (comm.rank() == 0) chains = nrg::discretize(gamma, sections[0].at("param"), z, longest, verbose ? "" : (fs::path(work) / "discretize.log").string());
    });
    broadcast(chains, comm);

    // Steps 6–8: each rank initializes and solves its twists
    auto const template_dir = fs::absolute(constr_params.get_model_dir()).string();
    char const *workdir_env = std::getenv("NRG_WORKDIR");
    collectively(comm, [&] {
      for (long i : mpi::chunk(range(long(z.size())), comm)) {
        auto const dir = fs::path(work) / std::to_string(i + 1);
        fs::create_directories(dir);
        std::ofstream(dir / "param") << section_text(sections[i]);
        auto const scratch = workdir_env ? fs::absolute(workdir_env).string() : dir.string();
        nrg::instantiate(dir.string(), template_dir, nmax[i], chains[i], verbose ? "" : (dir / "instantiate.log").string());
        nrg::solve(dir.string(), scratch, verbose ? "" : (dir / "solve.log").string());
      }
    });
    sync();

    // Raw results: every rank reads all twists (the working directory is shared)
    for (std::size_t k = 0; k < z.size(); ++k) {
      r.twists.push_back(read_twist((fs::path(work) / std::to_string(k + 1)).string(), z[k], spectra, options.fdm, options.energies > 0));
      r.twists.back().chain = chains[k];
    }
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
