#!/usr/bin/env python3

# The DMFT loop from Python with the NRG solver: the Hubbard model on the Bethe lattice at half filling, a few
# iterations; an interrupted run resumed from its checkpoint gives the iterations of an uninterrupted one

import shutil
import unittest

import numpy as np

from nrgljubljana_interface import Solver
from nrgljubljana_interface.dmft import BetheLattice, LinearMixer, Convergence, DMFTLoop, StepResultT

U = 1.0

SOLVE_KW = dict(T=1e-3, Lambda=4.0, Nz=2, Tmin=1e-4, keep=50, keepenergy=6.0, model_parameters={"U1": U})


class test_dmft_loop(unittest.TestCase):

    def make_loop(self, run_dir, alpha=0.5):
        S = Solver(model="SIAM", symtype="QS", mesh_max=4.0, mesh_min=1e-3, mesh_ratio=1.1)
        S.set_nrg_params(bins=50)
        loop = DMFTLoop(S, BetheLattice(1.0), LinearMixer(alpha), Convergence(eps=1e-8, max_iter=10), run_dir=run_dir, **SOLVE_KW)
        steps = []
        loop.set_callback(lambda step: steps.append(step))
        return S, loop, steps

    def Sigma0(self, S):
        Sigma = S.Sigma_w.copy()
        for name, s in Sigma:
            s.data[:] = U / 2  # the Hartree self-energy at half filling
        return Sigma

    def test_resume(self):
        for d in ("run_full", "run_resumed"):
            shutil.rmtree(d, ignore_errors=True)

        S, full, full_steps = self.make_loop("run_full")
        full.initialize(self.Sigma0(S), U / 2)
        self.assertEqual(full.run(3), "continue")
        self.assertEqual([s.iteration for s in full_steps], [1, 2, 3])
        self.assertIsInstance(full_steps[0], StepResultT)
        for s in full_steps:
            self.assertAlmostEqual(s.filling, 1.0, delta=1e-2)
            self.assertEqual(s.mu_next, U / 2)  # no μ control: μ stays fixed

        S1, first, _ = self.make_loop("run_resumed")
        first.initialize(U / 2, U / 2)  # a constant Σ0: the same as self.Sigma0
        first.run(2)

        S2, resumed, resumed_steps = self.make_loop("run_resumed")
        resumed.resume("run_resumed")
        self.assertEqual(resumed.iteration(), 2)
        resumed.run(1)
        self.assertEqual([s.iteration for s in resumed_steps], [3])
        np.testing.assert_allclose(resumed.Delta()['imp'].data, full.Delta()['imp'].data, atol=1e-8)
        self.assertEqual(resumed.mu(), full.mu())

        # another mixer: the checkpoint does not match
        S3, other, _ = self.make_loop("run_resumed", alpha=0.3)
        with self.assertRaisesRegex(Exception, "the mixer"):
            other.resume("run_resumed")

    def test_levels_are_set_by_the_loop(self):
        S = Solver(model="SIAM", symtype="QS", mesh_max=4.0, mesh_min=1e-3, mesh_ratio=1.1)
        kw = dict(SOLVE_KW, model_parameters={"U1": U, "eps1": -U / 2})
        with self.assertRaisesRegex(Exception, "the model parameter eps1 is an impurity level"):
            DMFTLoop(S, BetheLattice(1.0), LinearMixer(0.5), Convergence(), **kw)


if __name__ == '__main__':
    unittest.main()
