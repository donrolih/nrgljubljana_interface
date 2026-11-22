#!/usr/bin/env python3

# The in-process route (pipeline=library) against the script route, on the inputs of 1_SIAM_QS: SIAM/QS copied as a
# library template (its scripts removed), serial or on several ranks. The differences come from the route
# (mixchain with the integral method and Steffen interpolation, the workflow initializer, Nmax with bandrescale);
# they are printed for the record.

import os
import re
import shutil
import unittest

import numpy as np

from triqs.utility import mpi
from nrgljubljana_interface import Solver, SemiCircular

SCRIPTS = ("prepare", "discretize", "instantiate", "wilson", "Collect", "mmalog")


def library_templates(root):
    """A template root whose SIAM/QS runs in-process"""
    source = os.path.join(os.environ["NRGIF_TEMPLATE_DIR"], "SIAM", "QS")
    target = os.path.join(root, "SIAM", "QS")
    if mpi.is_master_node():
        shutil.rmtree(root, ignore_errors=True)
        shutil.copytree(source, target)
        for s in SCRIPTS:
            if os.path.exists(os.path.join(target, s)):
                os.remove(os.path.join(target, s))
        with open(os.path.join(target, "manifest")) as f:
            text = f.read()
        text = re.sub(r"^params=(.*)$", r"params=\1\npipeline=library", text, count=1, flags=re.M)
        with open(os.path.join(target, "manifest"), "w") as f:
            f.write(text)
    mpi.barrier()
    return root


def solve():
    S = Solver(model="SIAM", symtype="QS", mesh_max=1.0, mesh_min=1e-3, mesh_ratio=1.1)
    S.set_nrg_params(bins=50)
    S.Delta_w['imp'] << 0.1 * SemiCircular(1.0)
    S.solve(T=1e-3, Lambda=4.0, Nz=2, Tmin=1e-4, keep=50, keepenergy=6.0, model_parameters={"U1": 0.5, "eps1": -0.24})
    return S


class test_library_route(unittest.TestCase):

    def test_against_scripts(self):
        scripts = solve()

        shipped = os.environ["NRGIF_TEMPLATE_DIR"]
        os.environ["NRGIF_TEMPLATE_DIR"] = library_templates("library_templates_np%d" % mpi.size)
        try:
            library = solve()
        finally:
            os.environ["NRGIF_TEMPLATE_DIR"] = shipped

        r = library.last_result()
        self.assertEqual(len(r.twists), 2)
        for t in r.twists:
            self.assertEqual(len(t.chain.channels), 1)
            self.assertGreater(t.chain.channels[0].theta, 0)
            self.assertGreater(len(t.chain.channels[0].xi), 10)
            self.assertIn("density_interpolation=steffen", t.param)
        self.assertGreaterEqual(r.Gamma_projected_points, 0)
        self.assertEqual(len(scripts.last_result().twists[0].chain.channels), 0)  # the script route has no chain

        w = np.asarray(scripts.A_w['imp'].mesh.values())
        dA = np.max(np.abs(library.A_w['imp'].data[:, 0, 0] - scripts.A_w['imp'].data[:, 0, 0]))
        dSigma = np.max(np.abs(library.Sigma_w['imp'].data[:, 0, 0] - scripts.Sigma_w['imp'].data[:, 0, 0]))
        dn = abs(library.expv["n_d"] - scripts.expv["n_d"])
        if mpi.is_master_node():
            print(f"library - scripts: max |dA| = {dA:.3g} (max A = {np.max(scripts.A_w['imp'].data.real):.3g}), "
                  f"max |dSigma| = {dSigma:.3g}, |dn_d| = {dn:.3g}, Gamma raised at {r.Gamma_projected_points} points")
        # the same physics: the routes differ in the discretization only
        self.assertLess(dn, 1e-2)
        self.assertLess(dA, 0.05 * np.max(scripts.A_w['imp'].data.real))


if __name__ == '__main__':
    unittest.main()
