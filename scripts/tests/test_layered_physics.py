# Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
import unittest
from pathlib import Path
import subprocess

import numpy as np

from scripts.layered_physics import diagonal_observables, exact_stack


class LayeredPhysicsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root=Path(__file__).resolve().parents[2]
        subprocess.run(["make","-C","backend/build","layered-contract"],cwd=cls.root,check=True,capture_output=True)

    def test_cpp_and_independent_python_basis_observables_agree(self):
        output=subprocess.check_output([str(self.root/"backend/build/layered_contract.exe")],cwd=self.root,text=True)
        cpp=np.loadtxt(output.splitlines(),delimiter=",")
        _,py,_=diagonal_observables(2,2,2)
        for column,key in enumerate(("ice","rho4","plaquette","spin_q0_sq","spin_qpi_sq","dipole_q0_sq","dipole_qpi_sq"),start=1):
            np.testing.assert_allclose(cpp[:,column],py[key],atol=1e-15,rtol=0,err_msg=key)

    def test_layer_flip_preserves_local_quantities_and_exchanges_qz_channels(self):
        lx = ly = lz = 2
        states, observables, _ = diagonal_observables(lx, ly, lz)
        area = lx * ly
        flipped = states.copy()
        flipped[:, area:2*area] *= -1
        _, mapped, _ = diagonal_observables(lx, ly, lz, flipped)
        for key in ("ice", "rho4", "plaquette"):
            np.testing.assert_allclose(mapped[key], observables[key], atol=0, rtol=0)
        for q0, qpi in (("spin_q0_sq", "spin_qpi_sq"), ("dipole_q0_sq", "dipole_qpi_sq")):
            np.testing.assert_allclose(mapped[q0], observables[qpi], atol=1e-15, rtol=0)
            np.testing.assert_allclose(mapped[qpi], observables[q0], atol=1e-15, rtol=0)

    def test_exact_thermal_spectrum_and_qz_channels_obey_uniform_sign_mapping(self):
        common = dict(Lx=2, Ly=2, Lz=2, J0=1., J1=.3, J2=.025, hx=.2, beta=3.)
        plus, _, energies_plus = exact_stack({**common, "J3": .27})
        minus, _, energies_minus = exact_stack({**common, "J3": -.27})
        np.testing.assert_allclose(energies_plus, energies_minus, atol=2e-13, rtol=0)
        for key in ("E", "cv", "ice", "rho4", "plaquette"):
            self.assertAlmostEqual(plus[key], minus[key], places=12)
        for q0, qpi in (("spin_q0_sq", "spin_qpi_sq"), ("dipole_q0_sq", "dipole_qpi_sq")):
            self.assertAlmostEqual(plus[q0], minus[qpi], places=12)
            self.assertAlmostEqual(plus[qpi], minus[q0], places=12)

    def test_four_state_polarized_ice_proxy_has_zero_uniform_random_mean(self):
        _, values, _ = diagonal_observables(2, 2, 2)
        self.assertAlmostEqual(float(values["ice"].mean()), 3/8, places=14)
        self.assertAlmostEqual(float(values["rho4"].mean()), 0., places=14)
if __name__ == "__main__":
    unittest.main()
