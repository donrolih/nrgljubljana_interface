#!/usr/bin/env python3
# Solve options as keyword arguments, and the result of a solve: its fields, HDF5 round trips with and without the
# raw bins, a forward-only and a single-twist solve

import unittest

import numpy

from nrgljubljana_interface import Solver, SemiCircular, OptionsT, ResultT, broaden
from nrgljubljana_interface.solver import OPTIONS_FIELDS

from h5 import HDFArchive


def make_solver():
    S = Solver(model="SIAM", symtype="QS", mesh_max=1.0, mesh_min=1e-3, mesh_ratio=1.1)
    S.set_nrg_params(bins=50)
    S.Delta_w['imp'] << 0.1 * SemiCircular(1.0)
    return S


solve_params = dict(T=1e-3, Lambda=4.0, Nz=2, Tmin=1e-4, keep=50, keepenergy=6.0,
                    model_parameters={"U1": 0.5, "eps1": -0.24})


def max_diff(a, b):
    return max(numpy.max(numpy.abs(ga.data - gb.data)) for (_, ga), (_, gb) in zip(a, b))


class test_options_and_result(unittest.TestCase):

    def test_options_fields(self):
        names = {n for n in dir(OptionsT()) if not n.startswith('_')}
        self.assertEqual(names, set(OPTIONS_FIELDS))

    def test_default_solve(self):
        S = make_solver()
        self.assertIsNone(S.last_result())
        S.solve(**solve_params)
        r = S.last_result()
        self.assertIsInstance(r, ResultT)

        self.assertEqual(r.model, "SIAM")
        self.assertEqual(r.symtype, "QS")
        self.assertTrue(r.options.fdm)
        self.assertEqual(list(r.z), [0.5, 1.0])
        self.assertTrue(r.complete_z_average)
        self.assertEqual(r.correlator_type["A_d-A_d"], "d")
        self.assertIn("[correlators]", r.manifest)
        self.assertEqual(len(r.twists), 2)
        for t in r.twists:
            self.assertEqual(set(t.spectra), {"A_d-A_d", "self_d-A_d", "A_d-self_d", "self_d-self_d", "n_d-n_d", "sigma_d-sigma_d"})
            bins = t.spectra["A_d-A_d"]
            self.assertGreater(len(bins.frequency), 0)
            self.assertEqual(len(bins.frequency), len(bins.weight))
            self.assertEqual(t.custom.columns[0], "T")
            self.assertIn("n_d", t.customfdm)
            self.assertEqual(len(t.energy_flow), 0)
        self.assertIn("n_d", r.expv)
        self.assertEqual(max_diff(r.processed["G"], S.G_w), 0)
        self.assertEqual(max_diff(r.processed["Sigma"], S.Sigma_w), 0)

        # HDF5, with the bins: broadens to the same G
        with HDFArchive("10_options_and_result.out.h5", 'w') as ar:
            ar["r"] = r
            ar["without_bins"] = r.without_bins()
        with HDFArchive("10_options_and_result.out.h5", 'r') as ar:
            s = ar["r"]
            n = ar["without_bins"]
        self.assertIsInstance(s, ResultT)
        self.assertEqual(s.options.fdm, r.options.fdm)
        self.assertEqual(list(s.z), list(r.z))
        self.assertTrue(numpy.array_equal(s.twists[0].spectra["A_d-A_d"].weight, r.twists[0].spectra["A_d-A_d"].weight))
        self.assertEqual(max_diff(broaden(s, "G"), S.G_w), 0)

        # Without the bins: everything else, but no broadening
        self.assertEqual(len(n.twists), 2)
        self.assertEqual(len(n.twists[0].spectra), 0)
        self.assertEqual(max_diff(n.processed["G"], S.G_w), 0)
        with self.assertRaises(Exception):
            broaden(n, "G")

    def test_forward_only(self):
        S = make_solver()
        S.solve(**solve_params, fdm=False, energies=5)
        r = S.last_result()
        self.assertFalse(r.options.fdm)
        self.assertEqual(r.options.energies, 5)
        for t in r.twists:
            self.assertEqual(len(t.spectra), 0)
            self.assertEqual(len(t.customfdm), 0)
            self.assertGreater(len(t.energy_flow), 0)
            for step in t.energy_flow:
                self.assertTrue(0 < len(step) <= 5)
        level = r.twists[0].energy_flow[0][0]
        self.assertGreater(level.degeneracy, 0)
        self.assertEqual(len(r.processed), 0)
        self.assertIsNone(S.G_w)
        self.assertIsNone(S.Sigma_w)

    def test_single_twist(self):
        S = make_solver()
        S.solve(**solve_params, z=[0.5])
        r = S.last_result()
        self.assertEqual(list(r.z), [0.5])
        self.assertFalse(r.complete_z_average)
        self.assertEqual(len(r.twists), 1)
        self.assertIn("G", r.processed)

    def test_options_object(self):
        S = make_solver()
        from nrgljubljana_interface.solver_core import SolveParamsT
        S.solve(SolveParamsT(**solve_params), OptionsT(correlators=["G"], keep_bins=False))
        r = S.last_result()
        self.assertEqual(list(r.options.correlators), ["G"])
        self.assertEqual(set(r.processed), {"G"})
        self.assertEqual(len(r.twists[0].spectra), 0)
        self.assertIsNotNone(S.G_w)
        self.assertIsNone(S.F_l_w)
        self.assertIsNone(S.Sigma_w)


if __name__ == '__main__':
    unittest.main()
