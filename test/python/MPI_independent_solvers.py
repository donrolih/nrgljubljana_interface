#!/usr/bin/env python3
# Independent solvers with mpi=False: each rank solves its own parameters at the same time, with the same results as
# serial solves of the same parameters

import unittest

import numpy

from nrgljubljana_interface import Solver, SemiCircular
from triqs.utility import mpi

from h5 import HDFArchive


def solve(U):
    S = Solver(model="SIAM", symtype="QS", mesh_max=1.0, mesh_min=1e-3, mesh_ratio=1.1, mpi=False)
    S.set_nrg_params(bins=50)
    S.Delta_w['imp'] << 0.1 * SemiCircular(1.0)
    S.solve(T=1e-3, Lambda=4.0, Nz=2, Tmin=1e-4, keep=50, keepenergy=6.0,
            model_parameters={"U1": U, "eps1": -U / 2})
    return S


def U_of(rank):
    return 0.3 + 0.1 * rank


class test_independent_solvers(unittest.TestCase):

    def test_independent(self):
        rank = mpi.world.rank
        S = solve(U_of(rank))
        with HDFArchive(f"MPI_independent_solvers_rank{rank}.out.h5", 'w') as ar:
            ar["G_w"] = S.G_w
            ar["expv"] = S.expv
        mpi.barrier()

        if rank == 0:
            for other in range(mpi.world.size):
                with HDFArchive(f"MPI_independent_solvers_rank{other}.out.h5", 'r') as ar:
                    G, expv = ar["G_w"], ar["expv"]
                serial = solve(U_of(other))
                diff = max(numpy.max(numpy.abs(a.data - b.data)) for (_, a), (_, b) in zip(G, serial.G_w))
                self.assertEqual(diff, 0, f"rank {other}")
                self.assertEqual(dict(expv), dict(serial.expv), f"rank {other}")
        mpi.barrier()


if __name__ == '__main__':
    unittest.main()
