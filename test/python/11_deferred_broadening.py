#!/usr/bin/env python3
# Broadening the stored bins again: broaden(result, ...) and S.broaden(...)

import unittest

import numpy

from nrgljubljana_interface import Solver, SemiCircular, MeshReFreqLog, broaden


def max_diff(a, b):
    return max(numpy.max(numpy.abs(ga.data - gb.data)) for (_, ga), (_, gb) in zip(a, b))


class test_deferred_broadening(unittest.TestCase):

    def test_broaden(self):
        S = Solver(model="SIAM", symtype="QS", mesh_max=1.0, mesh_min=1e-3, mesh_ratio=1.1)
        S.set_nrg_params(bins=50)
        S.Delta_w['imp'] << 0.1 * SemiCircular(1.0)
        S.solve(T=1e-3, Lambda=4.0, Nz=2, Tmin=1e-4, keep=50, keepenergy=6.0,
                model_parameters={"U1": 0.5, "eps1": -0.24})
        r = S.last_result()
        G, F_l, Sigma = S.G_w.copy(), S.F_l_w.copy(), S.Sigma_w.copy()

        # Free function: the parameters of the solve give the solve's Green's functions; it changes nothing
        for name in ["G", "F_l", "F_r", "I", "chi_NN", "chi_SS"]:
            self.assertEqual(max_diff(broaden(r, name), r.processed[name]), 0, name)
        self.assertGreater(max_diff(broaden(r, "G", alpha=0.6), G), 1e-8)
        self.assertGreater(max_diff(broaden(r, "G", gamma=1.0), G), 1e-8)
        mesh = MeshReFreqLog(eps=1e-2, w_max=0.5, ratio=1.2)
        self.assertEqual(broaden(r, "G", mesh=mesh)['imp'].mesh, mesh)
        with self.assertRaises(Exception):
            broaden(r, "Sigma")
        self.assertEqual(max_diff(S.G_w, G), 0)

        # Method: updates the last result and G_w, A_w; the others keep their broadening
        g = S.broaden("G", alpha=0.6)
        self.assertGreater(max_diff(g, G), 1e-8)
        self.assertEqual(max_diff(S.G_w, g), 0)
        self.assertEqual(max_diff(S.last_result().processed["G"], g), 0)
        A = g.copy()
        for _, b in A:
            b.data[:] = -b.data.imag / numpy.pi
        self.assertEqual(max_diff(S.A_w, A), 0)
        self.assertEqual(max_diff(S.F_l_w, F_l), 0)
        self.assertEqual(max_diff(S.Sigma_w, Sigma), 0)

        # The parameters of the solve restore G
        S.broaden("G")
        self.assertEqual(max_diff(S.G_w, G), 0)


if __name__ == '__main__':
    unittest.main()
