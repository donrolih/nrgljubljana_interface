.. _aim:

Anderson impurity model (SIAM)
==============================

To illustrate how the NRG solver works in practice, we show the example of
one-orbital Anderson impurity embedded in a flat (Wilson) conduction bath. The
interacting part of the local Hamiltonian of this problem is simply:

.. math::

  \mathcal{H}_\mathrm{int} = U n_\uparrow n_\downarrow,

and the non-interacting Green's function is:

.. math::

    G^{-1}_{0,\sigma} (\omega) = \omega - \epsilon_f - V^2 g(\omega).

In this example, there is a Coulomb interaction :math:`U` on the impurity level,
which is at an energy :math:`\epsilon_f`. The bath Green's function is
:math:`g(\omega)`, and it has a flat density of states over the interval
:math:`[-1,1]`.  Finally, :math:`V` is the hybridization amplitude between the
impurity and the bath. Let us solve this problem with the NRG solver. Here is
the python :download:`script <aim.py>`:

.. literalinclude:: aim.py

Running this script takes a few minutes and generates an
HDF5 archive file called :file:`aim_solution.h5`. This file contains
the impurity spectral function, Green's function, self-energy, the expectation values of local variables, here
:math:`\langle n \rangle` and :math:`\langle n^2 \rangle`, and the result object of the solve. Let us plot
the spectral function:

.. plot:: guide/aim_plot.py
   :include-source:
   :scale: 100

As expected, the result shows a particle-hole symmetric impurity
spectral function.

Let us now go through the script in some more detail:

.. literalinclude:: aim.py
  :lines: 1-3

These lines import the classes to run the NRG solver, two functions of the DMFT layer
(:ref:`dmft`), and the class to store the results in a HDF5 archive file.

.. literalinclude:: aim.py
  :lines: 6-8

These lines define the impurity parameters, the (half)bandwidth :math:`D=1`,
and the physical temperature :math:`T`.

.. literalinclude:: aim.py
  :lines: 11

Here we construct the Solver object. The most important parameters here
are the model name and its symmetry type, as well as the mesh parameters which
define the logarithmic discretization mesh on which all quantities (input and
output) are defined. The models are described in template files in the directory
``templates/``. The program
will automatically determine the correct block structure of hybridisation
and Green's functions from the templates. This ensures that all quantities
within Solver are correctly initialized, and in particular that the block structure
of all Green's function objects is consistent.

.. literalinclude:: aim.py
  :lines: 14

This line defines, in addition to the physical temperature,
the key solve parameters which affect the quality of the results: the
discretization parameter :math:`\Lambda`, the number of the interleaved
discretization meshes :math:`N_z`, the scale ``Tmin`` which defines the length
of the Wilson chain (and which should be lower than the lowest physical temperature
scale in the problem, in this case lower than the Kondo temperature), and
the trucation settings which are determined by two parameters, the maximum number
of states kept ``keep`` and the truncation energy cutoff ``keepenergy``.
We also explicitly set the NRG energy window to the interval :math:`[-1:1]`, since
the default is to use an energy window with the same extent as the
logarithmic mesh, which in this example would be :math:`[-2:2]`, i.e., larger
than the support of the hybridisation function that we actually use; this
is achieved by setting the parameter ``bandrescale`` (back) to 1.
All other parameters have suitable default values for our purposes.

Finally, we have set the model parameters. These are contained in a separate Python
dictionary:

.. literalinclude:: aim.py
  :lines: 17-18

The hybridisation function is set by the following line:

.. literalinclude:: aim.py
  :lines: 21

We can now call the solver:

.. literalinclude:: aim.py
  :lines: 24

The solve options are given next to the solve parameters (see :ref:`options_results`). Here ``correlators``
selects the four correlators from which the self-energy is computed, :math:`\Sigma = \Sigma_H + I - F_l G^{-1} F_r`;
the template provides more (for example the charge and spin susceptibilities), which are then not computed.

The spectral bins of each correlator are kept in the result of the solve, so a correlator can be broadened again,
here with a wider broadening, without running NRG:

.. literalinclude:: aim.py
  :lines: 27

This replaces ``S.G_w`` and ``S.A_w``; ``S.Sigma_w`` is still that of the solve. The results are stored in an HDF5
archive, together with the result object of the solve (options, parameters, expectation values, the broadened
correlators). ``without_bins()`` leaves out the spectral bins, most of the size of a result; a result
stored with its bins can be broadened again later.

.. literalinclude:: aim.py
  :lines: 30-35

Finally, the occupancy computed from the Green's function, with the spin degeneracy of the block, is compared with
the expectation value ``S.expv["n_d"]`` from NRG, and the spectral weight of :math:`A(\omega)` is checked:

.. literalinclude:: aim.py
  :lines: 38-39

A version of this example can also be found in ``tutorial/1_AIM``
in the form of a Python script and a Jupyter (IPython) notebook.
