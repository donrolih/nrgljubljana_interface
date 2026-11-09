#!/usr/bin/env python3

# Hilbert transforms exposed by the package, on a refreq_log mesh

import unittest

import numpy as np

from nrgljubljana_interface import MeshReFreqLog, hilbert_transform_refreq, hilbert_transform_elementwise
from triqs.gfs import Gf


class test_hilbert_transform(unittest.TestCase):

    def setUp(self):
        self.mesh = MeshReFreqLog(eps=1e-4, w_max=1.0, ratio=1.01)
        self.z = 0.3 + 0.1j
        # Constant spectral matrix on [-1, 1]
        self.rho = np.array([[1.0, 0.5], [3.0, 2.0]])

    def test_scalar(self):
        A = Gf(mesh=self.mesh, target_shape=[])
        A.data[:] = 1.0
        # int_{-1}^{1} dw / (z - w) = log((z + 1) / (z - 1)); the gap (-1e-4, 1e-4) contributes O(1e-4)
        expected = np.log((self.z + 1) / (self.z - 1))
        self.assertAlmostEqual(hilbert_transform_refreq(A, self.z), expected, delta=1e-3)

    def test_elementwise(self):
        A = Gf(mesh=self.mesh, target_shape=[2, 2])
        A.data[:] = self.rho
        result = hilbert_transform_elementwise(A, self.z)
        self.assertEqual(result.shape, (2, 2))
        # Each element is the scalar transform of that element
        for i in range(2):
            for j in range(2):
                Aij = Gf(mesh=self.mesh, target_shape=[])
                Aij.data[:] = self.rho[i, j]
                self.assertAlmostEqual(result[i, j], hilbert_transform_refreq(Aij, self.z), places=12)


if __name__ == '__main__':
    unittest.main()
