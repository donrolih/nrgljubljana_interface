#!/usr/bin/env python3

# The DMFT functions and lattices from Python, on analytic input (no NRG run)

import unittest

import numpy as np

from h5 import HDFArchive
from triqs.gfs import Gf, BlockGf
from nrgljubljana_interface import MeshReFreqLog
from nrgljubljana_interface.dmft import (causal_projection, spectral_part, high_frequency_limit, rebuild_hermitian_part,
                                         hybridization_from_Gamma, Sigma_FG, Sigma_IFG, impurity_gf, occupancy,
                                         spectral_weight, spectral_distance, hilbert_transform, BetheLattice,
                                         TabulatedLattice, CallableLattice, DosT, LinearMixer, MuControl,
                                         Convergence, MuControlParamsT, ConvergenceParamsT, LoopParamsT, StepResultT,
                                         CausalityReportT, MuInputT, MixInputT)

D = 1.0


def bethe_g0(zeta):
    """G0 of the Bethe lattice of half-bandwidth D at ζ = ω + μ − Σ (upper half plane)"""
    return 2 / (zeta + np.sqrt(zeta - D) * np.sqrt(zeta + D))


def semicircle(n=4001):
    e = np.linspace(-D, D, n)
    return DosT(energy=list(e), density=list(2 / (np.pi * D**2) * np.sqrt(np.maximum(D**2 - e**2, 0))))


class test_dmft_functions(unittest.TestCase):

    def setUp(self):
        self.mesh = MeshReFreqLog(eps=1e-3, w_max=4.0, ratio=1.05)
        self.w = np.asarray(self.mesh.values())
        self.Sigma = self.block_gf(lambda w: 0.1 * w - 0.2j)

    def block_gf(self, f, names=('imp',)):
        blocks = []
        for _ in names:
            g = Gf(mesh=self.mesh, target_shape=[1, 1])
            g.data[:, 0, 0] = f(self.w)
            blocks.append(g)
        return BlockGf(name_list=list(names), block_list=blocks, make_copies=True)

    # --- Hilbert transform of a tabulated DOS

    def test_hilbert_transform(self):
        dos = semicircle()
        z = 0.3 + 0.1j
        self.assertAlmostEqual(hilbert_transform(dos, z), bethe_g0(z), delta=1e-4)
        # G(z̄) = conj G(z)
        self.assertAlmostEqual(hilbert_transform(dos, np.conj(z)), np.conj(hilbert_transform(dos, z)), places=14)
        # the first moment: G1 = (D²/4) G0²
        self.assertAlmostEqual(hilbert_transform(dos, z, 1), D**2 / 4 * bethe_g0(z)**2, delta=1e-4)
        zs = np.array([0.3 + 0.1j, -1.5 + 0.01j, 2.0j])
        np.testing.assert_allclose(hilbert_transform(dos, zs), bethe_g0(zs), atol=1e-4)

    def test_hilbert_transform_matrix(self):
        dos = semicircle()
        U = np.array([[1.0, 0.4j], [0.2, 1.0]])
        lam = np.array([0.3 + 0.1j, -0.5 + 0.2j])
        Z = U @ np.diag(lam) @ np.linalg.inv(U)
        expected = U @ np.diag([hilbert_transform(dos, l) for l in lam]) @ np.linalg.inv(U)
        np.testing.assert_allclose(hilbert_transform(dos, Z), expected, atol=1e-12)
        with self.assertRaisesRegex(Exception, "not finite"):
            hilbert_transform(dos, np.array([[np.nan, 0], [0, 1j]], dtype=complex))

    # --- lattices

    def test_lattice_classes(self):
        self.assertIn("bethe", BetheLattice(D).description())
        flat = DosT(energy=list(np.linspace(-D, D, 11)), density=[1 / (2 * D)] * 11)  # normalized exactly
        self.assertIn("tabulated shared", TabulatedLattice(flat).description())
        self.assertEqual(LinearMixer(0.5).alpha(), 0.5)

    def test_callable_lattice_equals_bethe(self):
        mu = 0.1

        def moments(Sigma, mu):
            G0, G1 = Sigma.copy(), Sigma.copy()
            for name, s in Sigma:
                zeta = self.w + mu - s.data[:, 0, 0]
                G0[name].data[:, 0, 0] = bethe_g0(zeta)
                G1[name].data[:, 0, 0] = D**2 / 4 * bethe_g0(zeta)**2
            return G0, G1

        callable_lattice = CallableLattice(moments, [np.zeros((1, 1), dtype=complex)], description="bethe in python")
        G_py, Delta_py = callable_lattice.self_consistency(self.Sigma, mu)
        G_cpp, Delta_cpp = BetheLattice(D).self_consistency(self.Sigma, mu)
        np.testing.assert_allclose(G_py['imp'].data, G_cpp['imp'].data, atol=1e-12)
        np.testing.assert_allclose(Delta_py['imp'].data, Delta_cpp['imp'].data, atol=1e-12)
        self.assertEqual(callable_lattice.description(), "bethe in python")

        # a moments function that fails raises its exception through the lattice
        def failing(Sigma, mu):
            raise ValueError("no moments here")

        with self.assertRaisesRegex(Exception, "no moments here"):
            CallableLattice(failing, [np.zeros((1, 1), dtype=complex)]).self_consistency(self.Sigma, mu)

    # --- causality, self-energy, measurements

    def test_causality(self):
        G, Delta = BetheLattice(D).self_consistency(self.Sigma, 0.0)
        projected, report = causal_projection(Delta, 1e-12)
        self.assertEqual(report.projected_points, 0)
        np.testing.assert_allclose(projected['imp'].data, Delta['imp'].data, atol=1e-14)
        S = spectral_part(Delta)
        self.assertTrue(np.all(S['imp'].data[:, 0, 0].real >= 0))
        self.assertEqual(len(high_frequency_limit(Delta)), 1)

        bad = Delta.copy()
        bad['imp'].data[10, 0, 0] += 0.5j  # a negative spectral part at one point
        _, report = causal_projection(bad, 1e-12)
        self.assertGreater(report.projected_points, 0)
        self.assertGreater(report.largest_violation, 0)

        # the Hermitian part from the spectral part by Kramers-Kronig, inside the band (the mesh truncates the tails)
        band = np.abs(self.w) < 1
        zero = [np.zeros((1, 1), dtype=complex)]
        rebuilt, _ = hybridization_from_Gamma(spectral_part(Delta), zero, 1e-12)
        self.assertLess(np.max(np.abs(rebuilt['imp'].data[band, 0, 0] - Delta['imp'].data[band, 0, 0])), 5e-2)
        rebuilt = rebuild_hermitian_part(Delta, zero)
        self.assertLess(np.max(np.abs(rebuilt['imp'].data[band, 0, 0] - Delta['imp'].data[band, 0, 0])), 5e-2)

    def test_self_energy(self):
        Delta = self.block_gf(lambda w: 0.2 * bethe_g0(w + 0.05j))
        h_imp = [np.array([[-0.3]], dtype=complex)]
        G = impurity_gf(Delta, self.Sigma, h_imp)
        expected = 1 / (self.w + 0.3 - Delta['imp'].data[:, 0, 0] - self.Sigma['imp'].data[:, 0, 0])
        np.testing.assert_allclose(G['imp'].data[:, 0, 0], expected, rtol=1e-12)

        F = G.copy()
        F['imp'].data[:] = self.Sigma['imp'].data * G['imp'].data
        np.testing.assert_allclose(Sigma_FG(G, F)['imp'].data, self.Sigma['imp'].data, rtol=1e-10)

        # Σ = Σ_H + I − F_l G⁻¹ F_r with F_l = F_r = 0 and I = Σ − Σ_H
        zero = G.copy()
        zero.zero()
        I = self.Sigma.copy()
        I['imp'].data[:] -= 0.25
        np.testing.assert_allclose(Sigma_IFG(G, zero, zero, I, [np.array([[0.25]], dtype=complex)])['imp'].data, self.Sigma['imp'].data, rtol=1e-12)

    def test_measurements(self):
        G = BetheLattice(D).local_gf(self.block_gf(lambda w: -1e-3j * np.ones_like(w)), 0.0)
        weight = spectral_weight(G)
        self.assertAlmostEqual(weight[0][0, 0].real, 1.0, delta=1e-2)
        n = occupancy(G, 1e-3, {'imp': 2})
        self.assertAlmostEqual(n.total, 1.0, delta=1e-2)
        self.assertEqual(spectral_distance(G, G).max, 0.0)

    # --- the loop controls one by one, as in a loop written in Python

    def test_modular_controls(self):
        lattice = BetheLattice(D)
        Sigma = self.block_gf(lambda w: -0.01j * np.ones_like(w))
        mu_control = MuControl(goal=0.8, maxdx=0.05, weight_tol=2e-2)
        adjusted = mu_control.adjust(MuInputT(Sigma=Sigma, mu=0.0, T=0.01, degeneracy={'imp': 2}), lattice)
        self.assertAlmostEqual(adjusted.filling, 1.0, delta=2e-2)  # half filling at μ = 0
        self.assertLess(adjusted.mu_next, 0.0)                     # towards the goal 0.8
        self.assertGreaterEqual(adjusted.mu_next, -0.05 - 1e-12)

        _, Delta = lattice.self_consistency(Sigma, 0.0)
        _, Delta_out = lattice.self_consistency(Sigma, adjusted.mu_next)
        Gamma, Gamma_out = spectral_part(Delta), spectral_part(Delta_out)
        mixed = LinearMixer(0.5).mix(MixInputT(iteration=1, Gamma=Gamma, Gamma_out=Gamma_out, mu=0.0, filling=adjusted.filling))
        np.testing.assert_allclose(mixed.Gamma_mix['imp'].data, 0.5 * (Gamma['imp'].data + Gamma_out['imp'].data), atol=1e-14)
        self.assertEqual(mixed.report.method, "linear")

    # --- parameters and reports in HDF5

    def test_h5_round_trip(self):
        mu_params = MuControlParamsT(goal=0.8, maxdx=0.02)
        conv_params = ConvergenceParamsT(eps=1e-5, window=4, eps_imp_loc=1e-3)
        loop_params = LoopParamsT(Gamma_floor=1e-5, run_dir="run", keep_history=False)
        report = CausalityReportT(projected_points=3, largest_violation=0.25)
        step = StepResultT(iteration=4, mu=0.1, mu_next=0.12, filling=0.8, kept_Delta=True)
        with HDFArchive("12_dmft_functions.out.h5", 'w') as ar:
            ar["mu_params"] = mu_params
            ar["conv_params"] = conv_params
            ar["loop_params"] = loop_params
            ar["report"] = report
            ar["step"] = step
        with HDFArchive("12_dmft_functions.out.h5", 'r') as ar:
            self.assertEqual(ar["mu_params"].goal, 0.8)
            self.assertEqual(ar["mu_params"].maxdx, 0.02)
            self.assertEqual(ar["conv_params"].window, 4)
            self.assertEqual(ar["conv_params"].eps_imp_loc, 1e-3)
            self.assertEqual(ar["loop_params"].Gamma_floor, 1e-5)
            self.assertEqual(ar["loop_params"].run_dir, "run")
            self.assertFalse(ar["loop_params"].keep_history)
            self.assertEqual(ar["report"].projected_points, 3)
            self.assertEqual(ar["report"].largest_violation, 0.25)
            self.assertEqual(ar["step"].iteration, 4)
            self.assertEqual(ar["step"].mu_next, 0.12)
            self.assertTrue(ar["step"].kept_Delta)

    # --- loop controls from keyword arguments, and C++ errors in Python

    def test_controls_and_errors(self):
        self.assertEqual(MuControl(goal=0.8, maxdx=0.01).params().maxdx, 0.01)
        self.assertEqual(Convergence(eps=1e-4, window=2).params().window, 2)
        with self.assertRaisesRegex(Exception, "α must satisfy 0 < α ≤ 1"):
            LinearMixer(0.0)
        with self.assertRaisesRegex(Exception, "window must be at least 1"):
            Convergence(window=0)


if __name__ == '__main__':
    unittest.main()
