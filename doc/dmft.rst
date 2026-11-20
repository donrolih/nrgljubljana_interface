.. _dmft:

The DMFT layer
==============

The package ``nrgljubljana_interface.dmft`` holds the pieces of a DMFT calculation with the NRG solver: functions
on retarded Green's functions, lattices, the loop controls, and a loop that combines them. The numerics are in C++;
the Python names are listed below. All Green's functions are ``BlockGf`` on the solver's logarithmic mesh, in the
block structure of the template; constant block matrices (on-site energies, Hartree self-energies, occupancies) are
lists of matrices in the block order.

Functions
---------

=================================================  ===============================================================
``spectral_part(F)``                               :math:`S = (F^\dagger - F)/2i`; for Δ this is Γ
``causal_projection(F, floor)``                    raises the eigenvalues of S to ``floor`` and corrects the
                                                   Hermitian part by Kramers-Kronig; returns ``(F, report)``
``rebuild_hermitian_part(F, F_inf)``               the Hermitian part from S and the high-frequency limit
``high_frequency_limit(F)``                        :math:`F_\infty` per block
``hybridization_from_Gamma(Gamma, h, floor)``      Δ from Γ and the on-site energy ⟨h⟩; returns ``(Delta, report)``
``Sigma_IFG(G, F_l, F_r, I, Sigma_H)``             :math:`\Sigma = \Sigma_H + I - F_l G^{-1} F_r`
``Sigma_FG(G, F_l)``                               :math:`\Sigma = F_l G^{-1}`
``impurity_gf(Delta, Sigma, h_imp)``               :math:`G = (\omega - h_{imp} - \Delta - \Sigma)^{-1}`
``occupancy(G, T, degeneracy)``                    occupancy matrices and their total, with the degeneracy of each
                                                   block (e.g. ``{"imp": 2}`` for spin)
``spectral_weight(G)``                             :math:`\int A\, d\omega` per block
``spectral_distance(a, b)``                        :math:`\int |A_a - A_b|\, d\omega` per element, and the largest
``hilbert_transform(dos, z, moment=0)``            :math:`\int \epsilon^n \rho(\epsilon)/(z - \epsilon)\, d\epsilon`
                                                   for a complex z, an array of z, or a square matrix z
=================================================  ===============================================================

Lattices
--------

A lattice gives, for a local self-energy and μ, the moments :math:`G_n = \frac{1}{N}\sum_k h_k^n (z - h_k)^{-1}`,
:math:`z = \omega + \mu - \Sigma`, and the on-site energy :math:`\langle h \rangle`. ``lattice.self_consistency(Sigma,
mu)`` returns the local Green's function and the hybridization :math:`\Delta = G_1 G_0^{-1}`.

- ``BetheLattice(half_bandwidth)``: closed forms.
- ``TabulatedLattice(dos)``: a tabulated DOS, ``DosT(energy=..., density=...)`` or
  ``TabulatedLattice.read_dos("DOS.dat")``, one for all orbitals or one per element; normalized to one.
- ``CallableLattice(moments, onsite_energy, local_gf=None, description="callable")``: a lattice written in Python;
  ``moments(Sigma, mu)`` returns ``(G0, G1)``. The description identifies the lattice when a run is resumed.

The loop
--------

``DMFTLoop(solver, lattice, mixer, convergence, mu_control=None, **kw)`` runs the iteration: the impurity levels
from ⟨h⟩ − μ, the solve, Σ made causal, μ from the μ control, the new Δ from the lattice, the convergence test, and
the mixing of Γ. The keyword arguments are those of ``Solver.solve`` (without the impurity levels) and
``Gamma_floor``, ``Sigma_floor``, ``run_dir``, ``keep_history``. :ref:`hubbard` is a complete example.

- ``LinearMixer(alpha)`` mixes Γ.
- ``MuControl(goal=..., maxdx=..., under=...)`` adjusts μ towards the filling ``goal`` with Σ held fixed; without a
  μ control, μ stays fixed.
- ``Convergence(eps=..., window=..., min_iter=..., max_iter=...)`` decides ``continue``, ``converged`` or ``stop``.

``loop.initialize(Sigma0, mu0)`` starts from a self-energy (a ``BlockGf`` or a constant), ``loop.step()`` runs one
iteration and returns its ``StepResultT``, ``loop.run(max_steps=0)`` runs until the status is ``converged`` or
``stop``. ``loop.set_callback(f)`` calls ``f(step)`` after each iteration. ``loop.mu()``, ``loop.Delta()``,
``loop.Sigma()``, ``loop.G_loc()`` and ``loop.last_step()`` give the current state.

With a run directory, rank 0 writes ``stats.dat`` (one line per iteration), ``checkpoint.h5`` (the latest
iteration, replaced atomically) and ``history.h5`` (all iterations). ``loop.resume(run_dir)`` continues from the
checkpoint; it checks that the parameters, the lattice, the mixer, the template and the mesh are those of the
checkpoint, and names each difference (``allow_changes=True`` continues with the new ones).
``loop.initialize_from(run_dir)`` starts a new run from the Δ and μ of another run.

Without the loop
----------------

The same steps can be written in Python, for full control of each one:

.. code-block:: python

   from nrgljubljana_interface import Solver
   from nrgljubljana_interface.dmft import (BetheLattice, MuControl, LinearMixer, MuInputT, MixInputT,
                                            causal_projection, spectral_part, hybridization_from_Gamma,
                                            spectral_distance)

   U, T = 4.0, 0.0025
   S = Solver(model="SIAM", symtype="QS", mesh_min=1e-7, mesh_max=40.0, mesh_ratio=1.01)
   lattice, mu_control, mixer = BetheLattice(1.0), MuControl(goal=0.8), LinearMixer(0.5)
   h = lattice.onsite_energy([("imp", 1)])

   mu = 0.51
   Sigma = S.Sigma_w.copy()
   Sigma << U / 2
   G_loc, Delta = lattice.self_consistency(Sigma, mu)
   for iteration in range(1, 51):
       S.Delta_w << Delta
       S.solve(T=T, Lambda=2.0, Nz=4, keep=2000, model_parameters={"U1": U, "eps1": h[0][0, 0].real - mu})
       Sigma, _ = causal_projection(S.Sigma_w, 1e-12)
       adjusted = mu_control.adjust(MuInputT(Sigma=Sigma, mu=mu, T=T, degeneracy={"imp": 2}), lattice)
       G_loc_next, Delta_out = lattice.self_consistency(Sigma, adjusted.mu_next)
       mixed = mixer.mix(MixInputT(iteration=iteration, Gamma=spectral_part(Delta), Gamma_out=spectral_part(Delta_out),
                                   mu=mu, filling=adjusted.filling))
       Delta, _ = hybridization_from_Gamma(mixed.Gamma_mix, h, 1e-6)
       distance = spectral_distance(G_loc_next, G_loc).max
       G_loc, mu = G_loc_next, adjusted.mu_next
       if distance < 1e-6:
           break
