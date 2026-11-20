# The semicircular density of states of the Bethe lattice, half-bandwidth D = 1, tabulated in DOS.dat: two columns,
# energy and density. In the band ε = D sin θ with θ uniform, which puts the points densely at the square-root band
# edges; outside the band ρ = 0 on a step of 0.001 up to |ε| = 1.3.
import numpy as np

D = 1.0
band = D * np.sin(np.linspace(-np.pi / 2, np.pi / 2, 2001))
outside = np.round(np.arange(1, 301) * 0.001, 12) + D
eps = np.concatenate([-outside[::-1], band, outside])
rho = 2 / (np.pi * D**2) * np.sqrt(np.clip(D**2 - eps**2, 0, None))
np.savetxt("DOS.dat", np.column_stack([eps, rho]), fmt="%.17g", delimiter="\t")
