#!/usr/bin/env python3

# The in-process route (pipeline=library) of SIAM/QS, serial or on several ranks: the chain coefficients of each
# twist, the rendered param, the projection of Gamma, and Nmax from Tmin in absolute units (with bandrescale).

import re
import unittest

import numpy as np

from nrgljubljana_interface import Solver, SemiCircular


def solve(bandrescale):
    S = Solver(model="SIAM", symtype="QS", mesh_max=1.0, mesh_min=1e-3, mesh_ratio=1.1)
    S.set_nrg_params(bins=50)
    S.Delta_w['imp'] << 0.1 * SemiCircular(1.0)
    S.solve(T=1e-3, Lambda=4.0, Nz=2, Tmin=1e-4, keep=50, keepenergy=6.0, bandrescale=bandrescale,
            model_parameters={"U1": 0.5, "eps1": -0.24})
    return S


def nmax(twist):
    return int(re.search(r"^Nmax=(\d+)$", twist.param, re.M).group(1))


class test_library_route(unittest.TestCase):

    def test_library_route(self):
        S = solve(1.0)
        r = S.last_result()
        self.assertEqual(len(r.twists), 2)
        for t in r.twists:
            self.assertEqual(len(t.chain.channels), 1)
            self.assertGreater(t.chain.channels[0].theta, 0)
            self.assertEqual(len(t.chain.channels[0].xi), len(t.chain.channels[0].zeta))
            self.assertGreaterEqual(len(t.chain.channels[0].xi), nmax(t))
            self.assertIn("density_interpolation=steffen", t.param)
            self.assertEqual(t.chain.Lambda, 4.0)
        self.assertEqual([t.chain.z for t in r.twists], [0.5, 1.0])
        self.assertGreater(r.Gamma_projected_points, 0)  # the zeros of the semicircle at its edges
        self.assertAlmostEqual(S.expv["n_d"], 1.0, delta=0.05)

        # Tmin is an absolute energy: with the band ten times wider in NRG units the chain reaches 10 Tmin in units
        # of the band edge, about 2 log_Lambda(10) = 3.3 sites more at Lambda = 4
        wide = solve(10.0).last_result()
        for narrow_twist, wide_twist in zip(r.twists, wide.twists):
            self.assertIn(nmax(wide_twist) - nmax(narrow_twist), (3, 4))


if __name__ == '__main__':
    unittest.main()
