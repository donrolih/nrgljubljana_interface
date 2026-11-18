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
#pragma once
#include "./container_set.hpp"
#include "./params.hpp"
#include "./types.hpp"
#include "./templates/manifest.hpp"
#include "./solver/options.hpp"
#include "./solver/result.hpp"
#include "./solver/process.hpp"

#include <iostream>
#include <fstream>
#include <map>
#include <string>

#include <triqs/gfs/hilbert_transform.hpp>
#include <triqs/utility/macros.hpp>

namespace nrgljubljana_interface {

  /**
   * @brief TRIQS interface to the NRGLjubljana numerical renormalization group impurity solver.
   *
   * @details Given a hybridization function, drives the external NRGLjubljana solver to compute
   * spectral functions, Green's functions, the self-energy, susceptibilities and
   * thermodynamic/expectation values for quantum impurity models.
   */
  class solver_core : public container_set {

    private:
    // The ranks that solve together: twists are split over them, and each holds the full result
    mpi::communicator comm;

    // Barrier on comm; no MPI call for a single rank
    void sync();

    // Fill the result containers (G_w, Sigma_w, expv, ...) from a processed result
    void fill_containers(result_t const &r);

    // The template's manifest, and the spectrum list (s, t, o, d, q) of each of its NRG correlators
    templates::manifest_t manifest;
    std::map<std::string, char> correlator_type;

    // The options of the current solve
    options_t options;

    // The result of the last solve
    std::optional<result_t> result;

    // Return reference to container_set
    container_set &result_set() { return static_cast<container_set &>(*this); }
    [[nodiscard]] container_set const &result_set() const { return static_cast<container_set const &>(*this); }

    public:
    /**
     * @brief Construct an NRGLjubljana_interface solver.
     *
     * @param cp Construction parameters.
     * @param mpi Solve on MPI_COMM_WORLD, splitting the twists over its ranks. If false, the solver makes no MPI
     * calls and solves all twists itself, so that independent solvers can run on different ranks.
     */
    explicit solver_core(constr_params_t cp, bool mpi = true);

    /**
     * @brief Construct a solver that works on the given communicator.
     *
     * @details The twists are split over its ranks. With a single rank (e.g. MPI_COMM_SELF) the solver makes no
     * collective MPI calls and solves all twists itself, so that independent solvers can run on different ranks.
     *
     * @param cp Construction parameters.
     * @param c The communicator.
     */
    C2PY_IGNORE solver_core(constr_params_t cp, mpi::communicator c);

    // Delete assignement operator because of const members
    solver_core(solver_core const &s)            = default;
    solver_core(solver_core &&s)                 = default;
    solver_core &operator=(solver_core const &s) = delete;
    solver_core &operator=(solver_core &&s)      = default;
    ~solver_core()                               = default;

    /**
     * @brief Solve the impurity problem.
     *
     * @details Runs the full NRGLjubljana pipeline for the assigned hybridization function and
     * populates the result containers.
     *
     * @param solve_params Parameters specific to the NRGLjubljana run.
     */
    void solve(solve_params_t const &solve_params);

    /**
     * @brief Solve the impurity problem with the given options.
     *
     * @details What NRG Ljubljana computes and what the result holds follow the options; solve(solve_params) is
     * this with the default options. A result container is set only for a correlator the solve broadened, except
     * that one the template does not provide is zero (if its structure is not empty).
     *
     * @param solve_params Parameters specific to the NRGLjubljana run.
     * @param options What to compute and return.
     */
    void solve(solve_params_t const &solve_params, options_t const &options);

    /**
     * @brief Broaden the bins of a correlator of the last result again, on the solver mesh, and update the solver.
     *
     * @details As broaden(last_result(), correlator, b), and the new Green's function is stored in the solver: it
     * replaces the correlator in last_result().processed and in its container (G_w for "G", F_l_w for "F_l", ...),
     * and the spectral density of that container is recomputed (A_w from G_w, B_l_w from F_l_w, ...). All other
     * containers are unchanged; in particular Sigma_w still comes from the broadening of the solve.
     *
     * @param correlator A correlator of the template.
     * @param b Broadening parameters; unset ones are those of the solve.
     * @return The Green's function of the correlator.
     */
    g_w_t broaden(std::string const &correlator, broadening_params_t const &b = {});

    /**
     * @brief Read and check the template's manifest, and set the block structures, the spectrum lists, the
     * model parameters and the polarization from it.
     *
     * @details Fails if the template has no manifest or the manifest does not match the template.
     */
    C2PY_IGNORE void read_manifest();

    /// The result of the last solve; std::nullopt before the first solve
    C2PY_IGNORE std::optional<result_t> const &last_result() const { return result; }

    /// The result of the last solve (a copy); None before the first solve
    C2PY_RENAME(last_result) std::optional<result_t> last_result_copy() const { return result; }

    /**
     * @brief Create a temporary working directory for a series of NRG runs.
     *
     * @param tempdir_ Base directory in which to create the temporary directory.
     * @return Path to the created temporary directory.
     */
    std::string create_tempdir(const std::string &tempdir_);

    /// Write \f$ \Gamma = -\mathrm{Im}\,\Delta(\omega) \f$ to a file.
    void write_gamma();

    /**
     * @brief Prepare the input files for an individual NRG calculation. Called from solve().
     *
     * @param z Discretization twist parameter.
     * @param taskdir Directory in which to prepare the calculation.
     */
    void instantiate(double z, const std::string &taskdir);

    /**
     * @brief Perform an individual NRG calculation. Called from ``solve()``.
     *
     * @param taskdir Directory containing the prepared calculation.
     */
    void solve_one(const std::string &taskdir);

    /**
     * @brief Adjust the advanced (low-level) NRG parameters.
     *
     * @param nrg_params Low-level NRG parameters.
     */
    void set_nrg_params(nrg_params_t const &nrg_params);

    /// Establish good defaults for the low-level NRG parameters.
    void set_params();

    /**
     * @brief Check that all required model parameters have been defined.
     *
     * @param sp Solve parameters to validate.
     */
    void check_model_params(const solve_params_t &sp);

    /**
     * @brief Produce the param file for a given value of the twist parameter \f$ z \f$.
     *
     * @param z Discretization twist parameter.
     */
    void generate_param_file(double z);

    /// Parameters used for the solver construction.
    constr_params_t constr_params;

    /// Low-level NRG parameters.
    nrg_params_t nrg_params;

    /// Parameters used for the most recent solve process.
    std::optional<solve_params_t> last_solve_params;

    /// If true, detailed output from NRGLjubljana and its tools is sent to stdout.
    bool verbose = false;

    /// Keep the temporary directories after the calculation.
    bool keep_temp_dir = false;

    /// The Green's function structure object.
    gf_struct_t gf_struct;

    /// The hybridization function structure object.
    gf_struct_t Delta_struct;

    /// The susceptibility structure object.
    gf_struct_t chi_struct;

    /// Logarithmic real-frequency mesh.
    refreq_log log_mesh;

    /// The hybridization function on the real-frequency axis.
    g_w_t Delta_w;

    /**
     * @brief Read expectation values from the NRG output files.
     *
     * @param Nz Number of discretization twists.
     */
    void readexpv(int Nz);

    /**
     * @brief Read thermodynamic variables (FDM algorithm) from the NRG output files.
     *
     * @param Nz Number of discretization twists.
     */
    void readtdfdm(int Nz);

    /// Read a scalar real-valued function name.dat
    // void readc(const std::string &name, std::optional<s_w_t> &s_w); // TO DO

    /// Suppress verbose output from the NRG solver.
    void be_quiet() { verbose = false; }

    /**
     * @brief Set the verbosity (see also ``be_quiet()``).
     *
     * @param v If true, enable verbose output.
     */
    void set_verbosity(bool v) { verbose = v; }

    /// HDF5 format tag for the solver object.
    static std::string hdf5_format() { return "NRGLJUBLJANA_INTERFACE_SolverCore"; }

    /// Write a solver object to an HDF5 file.
    friend void h5_write(h5::group h5group, std::string subgroup_name, solver_core const &s);

    /// Construct a solver object from an HDF5 file.
    C2PY_IGNORE
    static solver_core h5_read_construct(h5::group h5group, std::string subgroup_name);
  };

  /**
   * @brief Hilbert transform for refreq objects (scalar).
   *
   * @param gf Scalar real-frequency Green's function.
   * @param z Complex frequency at which to evaluate the transform.
   * @return Value of the Hilbert transform at \f$ z \f$.
   */
  std::complex<double> hilbert_transform_refreq(const c_w_cvt &gf, std::complex<double> z);

  /**
   * @brief Hilbert transform for refreq objects (matrix, elementwise).
   *
   * @param gf Matrix-valued real-frequency Green's function.
   * @param z Complex frequency at which to evaluate the transform.
   * @return Value of the elementwise Hilbert transform at \f$ z \f$.
   */
  matrix<std::complex<double>> hilbert_transform_elementwise(const m_w_cvt &gf, std::complex<double> z);

} // namespace nrgljubljana_interface
