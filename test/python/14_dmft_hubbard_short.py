#!/usr/bin/env python3

# The Hubbard model on the Bethe lattice with a filling goal, as doc/guide/hubbard.py at reduced resolution and with
# the closed forms of the Bethe lattice instead of the tabulated DOS: two iterations, serial or on several ranks

import os
import shutil
import unittest

from triqs.utility import mpi
from nrgljubljana_interface import Solver
from nrgljubljana_interface.dmft import BetheLattice, LinearMixer, MuControl, Convergence, DMFTLoop

U, goal, T = 2.0, 0.8, 0.01


class test_dmft_hubbard_short(unittest.TestCase):

    def test_two_iterations(self):
        run_dir = "run_hubbard_short_np%d" % mpi.size
        if mpi.is_master_node():
            shutil.rmtree(run_dir, ignore_errors=True)
        mpi.barrier()

        S = Solver(model="SIAM", symtype="QS", mesh_min=1e-5, mesh_max=10.0, mesh_ratio=1.03)
        S.set_nrg_params(bins=50)
        loop = DMFTLoop(S, BetheLattice(1.0), LinearMixer(0.5), Convergence(eps=1e-6, min_iter=5),
                        mu_control=MuControl(goal=goal, maxdx=0.05, weight_tol=1e-3), run_dir=run_dir,
                        T=T, Lambda=4.0, Nz=2, Tmin=1e-3, keep=100, keepenergy=6.0, model_parameters={"U1": U})
        steps = []
        loop.set_callback(lambda step: steps.append(step))
        loop.initialize(U / 2 - 1e-4j, U / 2)  # half filling
        self.assertEqual(loop.run(2), "continue")

        self.assertEqual([s.iteration for s in steps], [1, 2])
        for s in steps:
            self.assertIn(s.mu_report.status, ("within_tolerance", "bracketed", "limited"))
            self.assertLessEqual(abs(s.mu_next - s.mu), 0.05 + 1e-12)
        # from half filling towards the goal: μ decreases, the filling approaches the goal
        self.assertLess(steps[0].mu_next, steps[0].mu)
        self.assertLess(abs(steps[1].filling - goal), abs(steps[0].filling - goal))
        # every rank holds the same state
        self.assertEqual(mpi.bcast(loop.mu()), loop.mu())

        self.assertTrue(os.path.exists(os.path.join(run_dir, "checkpoint.h5")))
        with open(os.path.join(run_dir, "stats.dat")) as f:
            self.assertEqual(len([l for l in f if not l.startswith("#")]), 2)


if __name__ == '__main__':
    unittest.main()
