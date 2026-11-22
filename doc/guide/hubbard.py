# DMFT(NRG) for the Hubbard model on the Bethe lattice: a hole-doped metal at U = 4, n = 0.8, T = 0.0025 (D = 1)

import os
from h5 import HDFArchive
from triqs.utility import mpi
from nrgljubljana_interface import Solver
from nrgljubljana_interface.dmft import TabulatedLattice, LinearMixer, MuControl, Convergence, DMFTLoop

# Physical parameters
U = 4.0      # Hubbard repulsion
goal = 0.8   # band filling, both spins
T = 0.0025   # temperature

# The impurity solver: the model template and the logarithmic frequency mesh
S = Solver(model="SIAM", symtype="QS", mesh_min=1e-7, mesh_max=40.0, mesh_ratio=1.01)
S.set_nrg_params(Nmax=68, xmax=30, bins=500)

# The lattice: the semicircular DOS tabulated in DOS.dat (see dos.py)
lattice = TabulatedLattice(TabulatedLattice.read_dos("DOS.dat"))

# The loop controls: mixing of Γ = −Im Δ, the chemical potential, and the convergence test
mixer = LinearMixer(alpha=0.5)
mu_control = MuControl(goal=goal, maxdx=0.05, under=1.0)
convergence = Convergence(eps=1e-6, window=3, min_iter=10, max_iter=50)

# The loop; the keyword arguments are those of Solver.solve, without the impurity level eps1, which the loop sets
loop = DMFTLoop(S, lattice, mixer, convergence, mu_control=mu_control,
                Gamma_floor=1e-6, Sigma_floor=1e-12, run_dir="run",
                T=T, Lambda=2.0, Nz=4, Tmin=-1, keep=20000, keepmin=500, keepenergy=10.0, bandrescale=10.0,
                model_parameters={"U1": U})

# Start from the Hartree self-energy, or continue an interrupted run
if os.path.exists("run/checkpoint.h5"):
    loop.resume("run")
else:
    loop.initialize(U / 2 - 1e-4j, 0.51)

status = loop.run()

# Store the result (every rank holds it; rank 0 writes)
if mpi.is_master_node():
    with HDFArchive("solution.h5", "w") as ar:
        ar["status"] = status
        ar["mu"] = loop.mu()
        ar["Sigma"] = loop.Sigma()
        ar["G_loc"] = loop.G_loc()
        ar["Delta"] = loop.Delta()
        ar["last_step"] = loop.last_step()
