# Comparison against stored references with their recorded tolerance.
#
# Each reference X.ref.h5 has a file X.ref.json next to it. Its "precision" is the largest allowed
# max |a - b| of the Green's function data (the measure of h5diff); "reason" and "measured" record how
# the value was chosen. Tolerances are fixed when a reference is recorded, not adjusted to make a test pass.

import json
import os

from triqs.utility.h5diff import h5diff


def reference_tolerance(ref):
    """The recorded tolerance of the reference file ref (X.ref.h5 -> X.ref.json)"""
    root, ext = os.path.splitext(ref)
    if ext != '.h5':
        raise ValueError(f'Not an HDF5 reference: {ref}')
    with open(root + '.json') as f:
        return json.load(f)


def compare_to_reference(out, ref):
    """Compare the output file out with the reference file ref, at the reference's recorded precision"""
    h5diff(out, ref, precision=reference_tolerance(ref)['precision'])
