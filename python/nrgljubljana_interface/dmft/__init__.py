"""
The DMFT layer of nrgljubljana_interface.

The compiled module ``_dmft`` holds the numerics; this package gives its classes their Python names and adds
keyword-argument constructors for :class:`MuControl`, :class:`Convergence` and :class:`DMFTLoop`.

- Causality: ``spectral_part``, ``causal_projection``, ``rebuild_hermitian_part``, ``high_frequency_limit``,
  ``hybridization_from_Gamma``
- Self-energy: ``Sigma_IFG``, ``Sigma_FG``, ``impurity_gf``
- Measurements: ``occupancy``, ``spectral_weight``, ``spectral_distance``
- Lattices: :class:`BetheLattice`, :class:`TabulatedLattice`, :class:`CallableLattice`, and ``hilbert_transform`` of a
  tabulated DOS
- Loop controls: :class:`LinearMixer`, :class:`MuControl`, :class:`Convergence`
- The loop: :class:`DMFTLoop`

Block matrices (on-site energy, Hartree self-energy, occupancy) are lists of matrices in the block order of the
Green's functions.
"""
from .. import solver_core as _solver_core  # registers the solver types that _dmft uses
from ..solver import OPTIONS_FIELDS
from ..solver_core import SolveParamsT, OptionsT

from ._dmft import spectral_part, causal_projection, rebuild_hermitian_part, high_frequency_limit, hybridization_from_Gamma
from ._dmft import Sigma_IFG, Sigma_FG, impurity_gf
from ._dmft import occupancy, spectral_weight, spectral_distance, hilbert_transform
from ._dmft import CausalityReportT, OccupancyT, DistanceT, DosT
from ._dmft import MixInputT, MixReportT, MixResultT, MuControlParamsT, MuInputT, MuReportT, MuResultT
from ._dmft import ConvergenceParamsT, ConvergenceInputT, ConvergenceResultT, LoopParamsT, StepResultT
from ._dmft import BetheLatticeT as BetheLattice, TabulatedLatticeT as TabulatedLattice, CallableLatticeT as CallableLattice
from ._dmft import LinearMixer
from ._dmft import MuControlT as _MuControlT, ConvergenceT as _ConvergenceT, PyDmftLoop as _DmftLoop

# Fields of LoopParamsT other than the solve parameters and the options
LOOP_FIELDS = ('Gamma_floor', 'Sigma_floor', 'run_dir', 'keep_history')


class MuControl(_MuControlT):
    """
    The chemical potential between iterations: fixed without a goal, otherwise adjusted towards the goal with Σ held
    fixed. Takes a MuControlParamsT positionally, or its fields as keyword arguments (goal, maxdx, under, root_tol,
    mu_tol, max_evaluations, weight_tol, epsabs, epsrel).
    """

    def __init__(self, *args, **params_kw):
        _MuControlT.__init__(self, args[0] if args else MuControlParamsT(**params_kw))


class Convergence(_ConvergenceT):
    """
    The convergence test of the loop. Takes a ConvergenceParamsT positionally, or its fields as keyword arguments
    (eps, window, min_iter, max_iter, filling_tol, mu_tol, eps_imp_loc).
    """

    def __init__(self, *args, **params_kw):
        _ConvergenceT.__init__(self, args[0] if args else ConvergenceParamsT(**params_kw))


class DMFTLoop(_DmftLoop):
    """
    The DMFT loop around a Solver.

    Parameters
    ----------
    solver : Solver
        The impurity solver; its template must give the impurity levels in [levels].
    lattice : BetheLattice, TabulatedLattice or CallableLattice
        Copied into the loop.
    mixer : LinearMixer
        Copied into the loop, which holds the mixer's state.
    convergence : Convergence
    args : LoopParamsT, optional
        Loop parameters passed positionally, instead of the keyword arguments.
    mu_control : MuControl, optional
        Without one, μ stays fixed.
    params_kw : dict {'param':value}
        The keyword arguments of Solver.solve (solve parameters and options; no impurity level, which the loop sets),
        and Gamma_floor, Sigma_floor, run_dir, keep_history.

    The loop keeps references to the solver and the convergence test, which it uses in each step.
    """

    def __init__(self, solver, lattice, mixer, convergence, *args, mu_control=None, **params_kw):
        if args:
            params = args[0]
        else:
            loop_kw = {k: params_kw.pop(k) for k in LOOP_FIELDS if k in params_kw}
            options_kw = {k: params_kw.pop(k) for k in OPTIONS_FIELDS if k in params_kw}
            params = LoopParamsT(solve_params=SolveParamsT(**params_kw), options=OptionsT(**options_kw), **loop_kw)
        _DmftLoop.__init__(self, solver, lattice, mixer, convergence, params, mu_control)
        self._parts = (solver, convergence)


__all__ = ['spectral_part', 'causal_projection', 'rebuild_hermitian_part', 'high_frequency_limit', 'hybridization_from_Gamma',
           'Sigma_IFG', 'Sigma_FG', 'impurity_gf', 'occupancy', 'spectral_weight', 'spectral_distance', 'hilbert_transform',
           'BetheLattice', 'TabulatedLattice', 'CallableLattice', 'DosT', 'LinearMixer', 'MuControl',
           'Convergence', 'DMFTLoop', 'CausalityReportT', 'OccupancyT', 'DistanceT', 'MixInputT', 'MixReportT', 'MixResultT',
           'MuControlParamsT', 'MuInputT', 'MuReportT', 'MuResultT', 'ConvergenceParamsT', 'ConvergenceInputT',
           'ConvergenceResultT', 'LoopParamsT', 'StepResultT']
