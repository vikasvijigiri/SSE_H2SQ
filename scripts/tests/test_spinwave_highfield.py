# Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
import unittest

import numpy as np

from scripts.spinwave_highfield import (
    diagonal_energy,
    exact_odd_gap,
    expanded_pair_energy,
    harmonic_modes,
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


if __name__ == "__main__":
    unittest.main()
