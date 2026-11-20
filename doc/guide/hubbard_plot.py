import numpy as np
import matplotlib.pyplot as plt

from h5 import HDFArchive

with HDFArchive('solution.h5', 'r') as ar:
    print("status:", ar['status'], " mu =", ar['mu'], " filling =", ar['last_step'].filling)

    # Spectral function of the local lattice Green's function
    G_loc = ar['G_loc']['imp']
    w = np.asarray(G_loc.mesh.values())
    plt.plot(w, -G_loc.data[:, 0, 0].imag / np.pi)
    plt.xlim(-5, 5)
    plt.xlabel(r'$\omega$')
    plt.ylabel(r'$A(\omega)$')
    plt.show()
