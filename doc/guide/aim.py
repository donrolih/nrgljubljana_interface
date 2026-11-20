from nrgljubljana_interface import Solver, Flat
from nrgljubljana_interface.dmft import occupancy, spectral_weight
from h5 import HDFArchive

# Parameters
D, V, U = 1.0, 0.25, 1.0
e_f = -U/2.0 # particle-hole symmetric case
T = 1e-5

# Set up the Solver
S = Solver(model = "SIAM", symtype = "QS", mesh_max = 2.0, mesh_min = 1e-5, mesh_ratio = 1.01)

# Solve Parameters
sp = { "T": T, "Lambda": 2.0, "Nz": 4, "Tmin": 1e-6, "keep": 2000, "keepenergy": 10.0, "bandrescale": 1.0 }

# Model Parameters
mp = { "U1": U, "eps1": e_f }
sp["model_parameters"] = mp

# Initialize hybridization function
S.Delta_w['imp'] << V**2 * Flat(D)

# Solve the impurity model, computing the correlators that give the self-energy
S.solve(**sp, correlators = ["G", "F_l", "F_r", "I"])

# Broaden the spectral bins of G again, with another broadening width
S.broaden("G", alpha = 0.5)

# Store the Result, including the solve's result object without its spectral bins
with HDFArchive("aim_solution.h5", 'w') as arch:
    arch["A_w"] = S.A_w
    arch["G_w"] = S.G_w
    arch["Sigma_w"] = S.Sigma_w
    arch["expv"] = S.expv
    arch["result"] = S.last_result().without_bins()

# The occupancy from the Green's function (both spins) and its spectral weight, against the NRG expectation value
n = occupancy(S.G_w, T, {"imp": 2})
print("<n>=", S.expv["n_d"], " from G:", n.total, " weight of A:", spectral_weight(S.G_w)[0][0, 0].real)
