# Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
import unittest

import numpy as np

from scripts.spinwave_highfield import (
    bloch_adjacency,
    bloch_bond_shells,
    diagonal_energy,
    exact_odd_gap,
    expanded_pair_energy,
    high_symmetry_path,
    harmonic_modes,
    pair_couplings,
)


class HighFieldSpinWaveTests(unittest.TestCase):
    def test_pair_expansion_matches_hamiltonian_up_to_configuration_independent_shift(self):
        direct = diagonal_energy(2, 2, 2, 1., .3, .025, -.2)
        expanded = expanded_pair_energy(2, 2, 2, 1., .3, .025, -.2)
        residual = direct - expanded
        self.assertLess(float(np.ptp(residual)), 1e-12)

    def test_noninteracting_polarized_gap_is_two_k(self):
        for field in (0.3, 1., 5.):
            self.assertAlmostEqual(float(harmonic_modes(2, 2, 2, 0., 0., 0., field).min()), 2*field, places=12)
            self.assertAlmostEqual(exact_odd_gap(2, 2, 2, 0., 0., 0., 0., field), 2*field, places=12)

    def test_even_stack_j3_sign_is_spectral_mapping_in_ed_and_harmonic_theory(self):
        plus = harmonic_modes(2, 2, 2, .3, .025, .2, 4.)
        minus = harmonic_modes(2, 2, 2, .3, .025, -.2, 4.)
        np.testing.assert_allclose(plus, minus, atol=1e-12, rtol=0.)
        gap_plus = exact_odd_gap(2, 2, 2, 1., .3, .025, .2, 4.)
        gap_minus = exact_odd_gap(2, 2, 2, 1., .3, .025, -.2, 4.)
        self.assertAlmostEqual(gap_plus, gap_minus, places=11)

    def test_harmonic_error_decreases_in_the_high_field_limit(self):
        errors = []
        for field in (2., 4., 8., 16.):
            hp = float(harmonic_modes(2, 2, 2, .3, .025, .3, field).min())
            ed = exact_odd_gap(2, 2, 2, 1., .3, .025, .3, field)
            errors.append(abs(hp-ed)/ed)
        self.assertTrue(all(a > b for a, b in zip(errors, errors[1:])), errors)

    def test_bloch_matrix_reproduces_finite_real_space_pair_spectrum(self):
        lx = ly = lz = 4
        nsites = lx * ly * lz
        # Build the finite real-space adjacency independently from the Bloch transform.
        edges = pair_couplings(lx, ly, lz, .3, .025, .2)
        adjacency = np.zeros((nsites, nsites))
        for (i, j), coupling in edges.items():
            adjacency[i, j] += coupling
            adjacency[j, i] += coupling
        finite = np.linalg.eigvalsh(adjacency)

        # Fold the physical 4x4 momentum grid into the checkerboard primitive BZ.
        qpoints = set()
        for mx in range(lx):
            for my in range(ly):
                kx, ky = 2*np.pi*mx/lx, 2*np.pi*my/ly
                q1, q2 = (kx+ky) % (2*np.pi), (kx-ky) % (2*np.pi)
                qpoints.add((round(q1, 12), round(q2, 12)))
        shells = bloch_bond_shells(.3, .025, .2)
        bloch = []
        for q1, q2 in qpoints:
            for mz in range(lz):
                kz = 2*np.pi*mz/lz
                bloch.extend(np.linalg.eigvalsh(bloch_adjacency(q1, q2, kz, shells)))
        np.testing.assert_allclose(np.sort(bloch), finite, atol=1e-9, rtol=0.)

    def test_bloch_dispersion_obeys_even_stack_j3_sign_shift(self):
        shells_positive = bloch_bond_shells(.3, .025, .2)
        shells_negative = bloch_bond_shells(.3, .025, -.2)
        for q1, q2, kz in ((0., 0., 0.), (np.pi, 0., .37), (np.pi, np.pi, np.pi)):
            plus = np.linalg.eigvalsh(bloch_adjacency(q1, q2, kz + np.pi, shells_positive))
            minus = np.linalg.eigvalsh(bloch_adjacency(q1, q2, kz, shells_negative))
            np.testing.assert_allclose(plus, minus, atol=1e-12, rtol=0.)

    def test_high_symmetry_path_has_consistent_nodes_and_distances(self):
        qpoints, distance, ticks, labels = high_symmetry_path(points_per_segment=8)
        self.assertEqual(labels, [r"$\Gamma$", "X", "M", r"$\Gamma$", "Z", "R", "A", "Z"])
        self.assertEqual(len(ticks), len(labels))
        self.assertTrue(np.all(np.diff(ticks) > 0.))
        np.testing.assert_allclose(distance[0], ticks[0], atol=0.)
        np.testing.assert_allclose(distance[-1], ticks[-1], atol=1e-12)
        np.testing.assert_allclose(qpoints[0], np.zeros(3), atol=0.)
        np.testing.assert_allclose(qpoints[-1], [0., 0., np.pi], atol=0.)


if __name__ == "__main__":
    unittest.main()
