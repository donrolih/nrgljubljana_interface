###############################################################################
#
# nrgljubljana_interface: A TRIQS based impurity solver
#
# Copyright (c) 2019 The Simons foundation
#   authors: Nils Wentzell
#
# nrgljubljana_interface is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# nrgljubljana_interface is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# nrgljubljana_interface. If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################
"""
High-level interface to the NRGLjubljana impurity solver.

Defines the :class:`Solver` convenience subclass over the compiled
:class:`~nrgljubljana_interface.solver_core.SolverCore`, accepting the construction,
solve, and low-level NRG parameters as keyword arguments and importing
:mod:`triqs.utility.mpi` so that MPI is initialized before the solver runs.
"""
from .solver_core import SolverCore, ConstrParamsT, SolveParamsT, NrgParamsT, OptionsT, BroadeningParamsT
from .solver_core import broaden as _broaden

from triqs.gfs import *
from triqs.utility import mpi


# Fields of OptionsT, which Solver.solve takes as keyword arguments next to the solve parameters
OPTIONS_FIELDS = ('fdm', 'correlators', 'broaden', 'keep_bins', 'energies', 'z', 'keep_workdir')


def broaden(result, correlator, mesh=None, **params_kw):
    """
    Broaden the stored bins of a correlator of a result, without running NRG Ljubljana.

    Parameters
    ----------
    result : ResultT
        A solve result, e.g. ``S.last_result()`` or read from an HDF5 archive; it must hold the bins.
    correlator : str
        A correlator of the template (G, F_l, F_r, I, chi_NN, ...).
    mesh : MeshReFreqLog, optional
        The mesh; the solver mesh of the result if not given.
    params_kw : dict {'param':value}
        alpha, gamma, omega0_ratio; unset ones are those of the solve.

    Returns
    -------
    The Green's function of the correlator.
    """
    return _broaden(result, correlator, BroadeningParamsT(**params_kw), mesh)


# === The SolverCore Wrapper

class Solver(SolverCore):
    r"""
    NRGLjubljana impurity solver.

    Thin Python wrapper over the compiled
    :class:`~nrgljubljana_interface.solver_core.SolverCore`. Construct it with the model,
    symmetry type and frequency mesh, assign the hybridization function
    ``S.Delta_w['imp'] << ...``, then call :meth:`solve` and read the results
    (``S.A_w``, ``S.G_w``, ``S.Sigma_w``, ...).

    The :meth:`__init__`, :meth:`solve` and :meth:`set_nrg_params` methods each accept
    either keyword arguments (which are used to build the corresponding parameter object)
    or a single parameter object positionally, so a stored ``S.constr_params`` /
    ``S.last_solve_params`` / ``S.nrg_params`` can be passed straight back in.

    Parameters
    ----------
    args : ConstrParamsT, optional
        A construction-parameter object (e.g. a stored ``constr_params``) passed positionally,
        instead of the keyword arguments below.
    model : str, optional
        Impurity model to solve (selects a template directory), e.g. ``'SIAM'``. Default ``'SIAM'``.
    symtype : str, optional
        NRGLjubljana symmetry code, e.g. ``'QS'``, ``'QSZ'``, ``'ISO'``. Default ``'QS'``.
    mesh_max : float, optional
        Maximum frequency of the logarithmic mesh. Default ``10``.
    mesh_min : float, optional
        Minimum frequency of the logarithmic mesh. Default ``1e-4``.
    mesh_ratio : float, optional
        Common ratio of the geometric (logarithmic) frequency mesh. Default ``1.05``.
    mpi : bool, optional
        Solve on MPI_COMM_WORLD, splitting the twists over its ranks (default). If ``False``, the solver makes no
        MPI calls and solves all twists itself, so that independent solvers can run on different ranks.

    Notes
    -----
    See :class:`~nrgljubljana_interface.solver_core.ConstrParamsT` for the full list of
    construction parameters.
    """

    def __init__(self, *args, mpi=True, **params_kw):
        cp = args[0] if args else ConstrParamsT(**params_kw)
        SolverCore.__init__(self, cp, mpi)

    def solve(self, *args, **params_kw):
        """
        Solve the impurity problem.

        Parameters
        ----------
        args : SolveParamsT, optional
            A solve-parameter object (e.g. a stored ``last_solve_params``), optionally followed by an OptionsT.
        params_kw : dict {'param':value}
            Solve parameters (Lambda, Nz, T, model_parameters, ...; see SolveParamsT) and solve options (fdm,
            correlators, broaden, keep_bins, energies, z, keep_workdir; see OptionsT).
        """
        options_kw = {k: params_kw.pop(k) for k in OPTIONS_FIELDS if k in params_kw}
        sp = args[0] if args else SolveParamsT(**params_kw)
        options = args[1] if len(args) > 1 else OptionsT(**options_kw)
        return SolverCore.solve(self, sp, options)

    def broaden(self, correlator, *args, **params_kw):
        """
        Broaden the bins of a correlator of the last result again, on the solver mesh, and update the solver.

        The new Green's function is stored in the solver: it replaces the correlator in last_result().processed and
        in its container (G_w for "G", F_l_w for "F_l", ...), and the spectral density of that container is
        recomputed (A_w from G_w, B_l_w from F_l_w, ...). All other containers are unchanged; in particular Sigma_w
        still comes from the broadening of the solve.

        Parameters
        ----------
        correlator : str
            A correlator of the template (G, F_l, F_r, I, chi_NN, ...).
        args : BroadeningParamsT, optional
            A broadening-parameter object.
        params_kw : dict {'param':value}
            alpha, gamma, omega0_ratio; unset ones are those of the solve.

        Returns
        -------
        The Green's function of the correlator.
        """
        b = args[0] if args else BroadeningParamsT(**params_kw)
        return SolverCore.broaden(self, correlator, b)

    def set_nrg_params(self, *args, **params_kw):
        """
        Set the advanced (low-level) NRG parameters.

        Parameters
        ----------
        args : NrgParamsT, optional
            A low-level NRG-parameter object (e.g. a stored ``nrg_params``).
        params_kw : dict {'param':value}
            Low-level NRG parameters; see NrgParamsT.
        """
        np = args[0] if args else NrgParamsT(**params_kw)
        return SolverCore.set_nrg_params(self, np)
