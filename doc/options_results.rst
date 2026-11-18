.. _options_results:

Solve options and results
=========================

.. code-block:: python

   from nrgljubljana_interface import Solver, broaden
   from h5 import HDFArchive

   S = Solver(model="SIAM", symtype="QS", mesh_min=1e-4, mesh_max=2.0, mesh_ratio=1.05)
   S.Delta_w['imp'] << ...
   S.solve(T=1e-3, Lambda=2.0, Nz=4, model_parameters={"U1": 0.5, "eps1": -0.25},
           correlators=["G", "F_l", "F_r", "I"])          # solve options next to the solve parameters

   r = S.last_result()                                    # options, parameters, raw data of each twist, ...
   with HDFArchive("result.h5", 'w') as ar:
       ar["result"] = r                                   # or r.without_bins(): smaller, cannot be broadened again

   G = broaden(r, "G", alpha=0.5)                         # broaden the stored bins again
   S.broaden("G", alpha=0.5)                              # the same, and update S.G_w and S.A_w

- **Options** (:class:`OptionsT <nrgljubljana_interface.solver_core.OptionsT>`): ``fdm``, ``correlators``,
  ``broaden``, ``keep_bins``, ``energies``, ``z``, ``keep_workdir``. The defaults compute and broaden all correlators
  of the template, as before. A container (``S.G_w``, ...) is set only for a correlator that was broadened.
- **Result** (:class:`ResultT <nrgljubljana_interface.solver_core.ResultT>`): what was solved, the raw data of each
  twist (spectral bins, ``custom``, ``td``, FDM values, energy flow), their averages and the broadened correlators.
- **Broadening** (:class:`BroadeningParamsT <nrgljubljana_interface.solver_core.BroadeningParamsT>`): ``alpha``,
  ``gamma`` (unset: those of the solve) and ``omega0_ratio``; T is that of the solve. Only correlators of the
  template, not Sigma.
- **Independent solvers**: ``Solver(..., mpi=False)`` makes no MPI calls and solves all twists itself, so that
  different MPI ranks can solve different problems at the same time.
