.. _hubbard:

DMFT for the Hubbard model
==========================

We now present a Python script for DMFT calculations with the NRG as the impurity solver: a strongly interacting,
hole-doped Hubbard model on the Bethe lattice, :math:`U = 4`, filling :math:`n = 0.8` (both spins) and
:math:`T = 0.0025`, in units of the half-bandwidth :math:`D = 1`. The DMFT loop, the lattice and the loop controls
come from the DMFT layer, :ref:`dmft`. Here is the :download:`listing <hubbard.py>`:

.. literalinclude:: hubbard.py

Running this script takes about an hour on eight cores. The NRG solver produces big temporary files during the
calculation; they are stored in a working directory created in the current directory. To override this, set the
environment variable ``NRG_WORKDIR``, for example to a local scratch disk.

Let us go through the script. The impurity solver is constructed as for a single impurity problem; the mesh must
contain the whole spectrum, here up to :math:`|\omega| = 40`:

.. literalinclude:: hubbard.py
  :lines: 14-15

The length of the Wilson chain is set directly by ``Nmax`` (with ``Tmin=-1`` below, so that it is not derived from
``Tmin``). The lattice is given by its density of states, tabulated in :download:`DOS.dat <DOS.dat>`, which the
script :download:`dos.py <dos.py>` generates; the DOS must be normalized to one. ``BetheLattice`` gives the same
lattice through closed forms, and ``CallableLattice`` accepts a lattice written in Python.

.. literalinclude:: hubbard.py
  :lines: 18

The loop controls: linear mixing of :math:`\Gamma = -\mathrm{Im}\,\Delta`, the chemical potential adjusted towards
the filling with the self-energy held fixed (at most 0.05 per iteration), and the convergence test, which requires
three consecutive iterations in which the local spectral function changes by less than :math:`10^{-6}`, between the
10th and the 50th iteration:

.. literalinclude:: hubbard.py
  :lines: 21-23

The loop takes the solver, the lattice and the loop controls, the floors of the spectral parts of :math:`\Delta`
and :math:`\Sigma`, the run directory, and the keyword arguments of ``Solver.solve``. The impurity level ``eps1``
is not among them: the template maps it to the impurity (its ``[levels]`` section), and the loop sets it from the
on-site energy of the lattice and :math:`\mu` in each iteration.

.. literalinclude:: hubbard.py
  :lines: 26-29

The loop starts from the Hartree self-energy :math:`U/2` and :math:`\mu = 0.51`, or continues an interrupted run
from the checkpoint in its run directory:

.. literalinclude:: hubbard.py
  :lines: 32-37

During the run the loop writes three files in ``run/``: ``stats.dat``, one line per iteration (μ and μ for the next
iteration, the filling, the convergence measures, the mixing residual, the causal projections, the status),
``checkpoint.h5``, the latest iteration, and ``history.h5``, all iterations. ``run()`` returns ``"converged"``, or
``"stop"`` if the loop reached ``max_iter``. The script then stores the result:

.. literalinclude:: hubbard.py
  :lines: 40-46

Let us plot the spectral function of the local lattice Green's function:

.. plot:: guide/hubbard_plot.py
   :include-source:
   :scale: 100
