# Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
"""Independent physical limits, response identities and correlated-data checks."""
import os
from pathlib import Path
import shutil
import subprocess
import sys
import unittest
import tempfile
import threading
import math
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "backend"))
from app.services.physics_model import exact_reference, spectrum, geometry
from app.services.benchmark_statistics import autocorrelation, split_rhat


class ReferenceTests(unittest.TestCase):
    def params(self, **values):
        return dict(Lx=4, Ly=2, J0=1., J1=.3, J2=.025, hx=.2, beta=2., **values)

    def test_noninteracting_quantum_limit(self):
        p = dict(Lx=2, Ly=2, J0=0., J1=0., J2=0., hx=.37, beta=1.7)
        v = exact_reference(p)["values"]
        expected = np.tanh(p["beta"]*p["hx"])
        self.assertAlmostEqual(v["mx"], expected, places=12)
        self.assertAlmostEqual(v["E"], -p["hx"]*expected, places=12)
        self.assertAlmostEqual(v["chi_x"], p["beta"]*(1-expected**2), places=12)
        self.assertAlmostEqual(v["mz2"], 1/16, places=12)

    def test_transverse_response_is_energy_derivative(self):
        p = self.params(); h = 1e-4
        center = exact_reference(p)
        plus = exact_reference({**p, "hx": p["hx"]+h})["values"]
        minus = exact_reference({**p, "hx": p["hx"]-h})["values"]
        self.assertAlmostEqual((plus["mx"]-minus["mx"])/(2*h), center["values"]["chi_x"], places=6)
        self.assertLess(center["numerical_residual"], 1e-11)

    def test_classical_complete_degenerate_manifold(self):
        p = dict(Lx=4, Ly=4, J0=1., J1=.3, J2=0., hx=0., beta=10.)
        v = exact_reference(p)["values"]
        self.assertAlmostEqual(v["E"], -.8, places=8)
        self.assertLess(v["mz2"], 1e-9)
        self.assertAlmostEqual(v["ice"], 1, places=8)
        self.assertAlmostEqual(v["ct0"], v["ct4"], places=14)

    def test_neighbor_sign_and_directional_multiplicity(self):
        for lx, ly in [(2, 2), (2, 4), (4, 2), (4, 4)]:
            plaquettes, pairs, signs = geometry(lx, ly)
            self.assertEqual(len(pairs), lx*ly)
            self.assertTrue(all(signs[a]*signs[b] == -1 for a, b in pairs))
            self.assertTrue(all(len(set(p)) == 4 for p in plaquettes))

    def test_quantum_reference_never_silently_truncates(self):
        with self.assertRaisesRegex(ValueError, "10 spins"):
            spectrum(4, 4, 1., .3, 0., .2)


class StatisticsTests(unittest.TestCase):
    def test_ar1_known_integrated_autocorrelation(self):
        rng = np.random.default_rng(178)
        noise = rng.normal(size=200000); x = np.zeros(len(noise))
        for i in range(1, len(x)): x[i] = .8*x[i-1]+noise[i]
        value = autocorrelation(x)
        self.assertAlmostEqual(value["tau"], 4.5, delta=.5)
        self.assertAlmostEqual(value["ess"], len(x)/9, delta=3000)

    def test_constant_chains_do_not_fake_zero_autocorrelation(self):
        self.assertIsNone(autocorrelation(np.ones(1000))["tau"])
        self.assertIsNone(autocorrelation(np.ones(1000))["ess"])
        self.assertTrue(np.isinf(split_rhat([np.ones(1000), np.zeros(1000)])))

    def test_rhat_detects_between_chain_disagreement(self):
        rng = np.random.default_rng(42)
        chains = [rng.normal(size=2000) for _ in range(4)]
        self.assertLess(split_rhat(chains), 1.01)
        chains[-1] += 2
        self.assertGreater(split_rhat(chains), 1.1)

    def test_exact_reference_catches_a_shared_frozen_sector(self):
        from app.services.benchmark import compare_reference
        summary = {"total_seconds": 1, "observables": {"ct0": {"mean": .001, "error": .00004, "ess": 1000, "status": "checks passed"}}}
        compare_reference(summary, {"values": {"ct0": .5}})
        self.assertEqual(summary["observables"]["ct0"]["status"], "review")
        self.assertLess(summary["observables"]["ct0"]["pull"], -1000)

    def test_saved_comparison_preserves_original_provenance_and_windows_paths(self):
        from app.services import benchmark
        params = dict(Lx=2, Ly=2, beta=1., J0=0., J1=0., J2=0., hx=0.)
        def fixture(identity, method, path):
            return {"id": identity, "created": "2026-10-16", "status": "complete", "request": {**params, "chains": 2, "n_therm": 100, "n_measure": 1000, "seed": 1, "axes": {"beta": [1.]}, "algorithms": [method]}, "completed_runs": 2, "total_runs": 2, "points": [{"axis": "beta", "x": 1., "point_id": "0"*16, "params": params, "ed": {"values": {"E": 0.}}, "methods": {method: {"observables": {"E": {"mean": 0., "error": .1, "status": "checks passed"}}}}}], "provenance": {"model_version": "test", "source_sha256": {path: "reference-hash"}}, "definitions": {}, "citations": []}
        with tempfile.TemporaryDirectory(dir=ROOT / ".cache") as directory, patch.object(benchmark, "STORE", Path(directory)):
            benchmark.save_job(fixture("1"*32, "vertex", "backend\\app\\services\\physics_model.py"))
            benchmark.save_job(fixture("2"*32, "tempered", "backend/app/services/physics_model.py"))
            combined = benchmark.combine_jobs(["1"*32, "2"*32])
            self.assertEqual(combined["completed_runs"], 4)
            self.assertEqual(set(combined["points"][0]["methods"]), {"vertex", "tempered"})
            self.assertEqual(combined["source_experiments"][0]["request"]["algorithms"], ["vertex"])
            self.assertEqual(benchmark.get_job("1"*32)["request"]["algorithms"], ["vertex"])
            self.assertEqual(combined["points"][0]["methods"]["tempered"]["source_experiment"], "2"*32)


class KernelIntegrationTests(unittest.TestCase):
    def test_quantum_measurements_and_closed_worldlines_against_complete_ed(self):
        from app.services.benchmark import EXECUTABLE, BenchmarkRequest, run_chain
        from app.services.benchmark_statistics import summarize
        if not EXECUTABLE.exists(): self.skipTest("Build the benchmark engine first")
        request = BenchmarkRequest(beta=1., J2=.05, hx=.3, chains=2, n_therm=2000, n_measure=8000, axes={"beta": [1.]})
        p = {key: getattr(request, key) for key in ("Lx", "Ly", "J0", "J1", "J2", "hx", "beta")}
        with tempfile.TemporaryDirectory(dir=ROOT / ".cache") as directory:
            runs = [run_chain(p, "canonical", 8700+i, Path(directory)/str(i), request, threading.Event()) for i in range(2)]
            self.assertTrue(all(meta["configuration_checks"] >= 10 for _, meta in runs))
            self.assertTrue(all(meta["cutoff_saturations"] == 0 for _, meta in runs))
            summary = summarize([trace for trace, _ in runs], p)
        reference = exact_reference(p)["values"]
        for key in ("E", "E2", "mx", "ice", "stripe2", "ct0", "chi_x"):
            with self.subTest(observable=key):
                value = summary["observables"][key]
                self.assertLess(abs(value["mean"]-reference[key]), 8*value["error"]+5e-4)


@unittest.skipUnless(shutil.which("g++"), "C++ compiler unavailable")
class CrossLanguageHamiltonianTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.exe = ROOT / ".cache" / ("diagonal_contract.exe" if os.name == "nt" else "diagonal_contract")
        cls.exe.parent.mkdir(exist_ok=True)
        sources = ["backend/benchmarks/diagonal_contract.cpp", "backend/src/variables.cpp", "backend/src/updates.cpp", "backend/src/observables.cpp"]
        command = ["g++", "-O2", "-std=c++17", *sources, "-o", str(cls.exe)]
        if os.name == "nt": command += ["-static-libgcc", "-static-libstdc++", "-Wl,-Bstatic,--whole-archive", "-lwinpthread", "-Wl,--no-whole-archive,-Bdynamic"]
        subprocess.run(command, cwd=ROOT, check=True, capture_output=True)

    def test_every_classical_basis_energy_in_rectangular_and_square_cells(self):
        for lx, ly, j0, j1, j2 in [(2, 2, 1., .3, .025), (4, 2, 1., .3, .05), (2, 4, 1., -.3, -.05), (4, 4, 1., .3, .025), (4, 2, -1., -.2, .03)]:
            with self.subTest(lattice=(lx, ly), couplings=(j0, j1, j2)):
                output = subprocess.check_output([str(self.exe), str(lx), str(ly), str(j0), str(j1), str(j2)], text=True)
                actual = np.fromstring(output, sep='\n')
                expected = spectrum(lx, ly, j0, j1, j2, 0.)[0]
                np.testing.assert_allclose(actual, expected, atol=3e-13, rtol=0)

    def test_stacked_basis_energies_match_independent_periodic_3d_hamiltonian(self):
        # Includes Lz=2, where the periodic +z sum contains two parallel
        # layer links, and rectangular cells to expose indexing mistakes.
        for lx, ly, lz, j0, j1, j2, j3 in [
            (2, 2, 2, 1., .3, .07, -.19),
            (4, 2, 2, 1., -.2, .035, .11),
            (2, 4, 2, -.8, .15, -.04, .23),
        ]:
            with self.subTest(lattice=(lx, ly, lz)):
                output = subprocess.check_output([
                    str(self.exe), "stack", str(lx), str(ly), str(lz),
                    str(j0), str(j1), str(j2), str(j3)
                ], text=True)
                actual = np.fromstring(output, sep="\n")

                area, n = lx*ly, lx*ly*lz
                plaquettes, pairs, _ = geometry(lx, ly)
                states = 2*((np.arange(1 << n)[:, None] >> np.arange(n)) & 1)-1
                expected = np.zeros(len(states))
                for z in range(lz):
                    layer = states[:, z*area:(z+1)*area]
                    local = layer[:, plaquettes]
                    expected -= j0*np.prod(local, axis=2).sum(axis=1)
                    expected += j1*(local[:, :, 0]*local[:, :, 2] +
                                    local[:, :, 1]*local[:, :, 3]).sum(axis=1)
                    px, py = local @ (np.array([1, -1, -1, 1])/4), local @ (np.array([1, 1, -1, -1])/4)
                    expected += j2*sum(px[:, a]*px[:, b] + py[:, a]*py[:, b] for a, b in pairs)
                for z in range(lz):
                    expected += j3*np.sum(states[:, z*area:(z+1)*area] *
                                          states[:, ((z+1) % lz)*area:((z+1) % lz+1)*area], axis=1)
                np.testing.assert_allclose(actual, expected, atol=4e-13, rtol=0)

    def test_uniform_interlayer_signs_have_identical_small_stack_spectra(self):
        # Exact finite-system check of the layer-flip unitary, including all
        # in-plane interactions and the noncommuting transverse-field term.
        lx = ly = lz = 2
        area, n = lx*ly, lx*ly*lz
        beta_params = dict(J0=1., J1=.27, J2=.041, hx=.23, J3=.19)
        plaquettes, pairs, _ = geometry(lx, ly)
        basis = np.arange(1 << n)
        spins = 2*((basis[:, None] >> np.arange(n)) & 1)-1

        def hamiltonian(j3):
            diagonal = np.zeros(len(basis))
            for z in range(lz):
                local = spins[:, z*area:(z+1)*area][:, plaquettes]
                diagonal -= beta_params["J0"]*np.prod(local, axis=2).sum(axis=1)
                diagonal += beta_params["J1"]*(local[:, :, 0]*local[:, :, 2] +
                                                 local[:, :, 1]*local[:, :, 3]).sum(axis=1)
                px = local @ (np.array([1, -1, -1, 1])/4)
                py = local @ (np.array([1, 1, -1, -1])/4)
                diagonal += beta_params["J2"]*sum(
                    px[:, a]*px[:, b] + py[:, a]*py[:, b] for a, b in pairs
                )
            for z in range(lz):
                diagonal += j3*np.sum(spins[:, z*area:(z+1)*area] *
                                      spins[:, ((z+1) % lz)*area:((z+1) % lz+1)*area], axis=1)
            matrix = np.diag(diagonal)
            for site in range(n):
                matrix[basis, basis ^ (1 << site)] = -beta_params["hx"]
            return matrix

        h_plus, h_minus = hamiltonian(beta_params["J3"]), hamiltonian(-beta_params["J3"])
        odd_layer_mask = sum(1 << (area + site) for site in range(area))
        transformed_basis = basis ^ odd_layer_mask
        np.testing.assert_allclose(
            h_plus[np.ix_(transformed_basis, transformed_basis)], h_minus,
            atol=2e-15, rtol=0
        )
        np.testing.assert_allclose(
            np.linalg.eigvalsh(h_plus), np.linalg.eigvalsh(h_minus), atol=2e-13, rtol=0
        )

    def test_large_stacked_geometry_indices_and_checkerboard_periods(self):
        for lx, ly, lz in [(8, 6, 4), (12, 10, 8)]:
            with self.subTest(lattice=(lx, ly, lz)):
                output = subprocess.check_output(
                    [str(self.exe), "geometry", str(lx), str(ly), str(lz)], text=True
                )
                ns, black, white, pairs, _bonds, vertical = map(int, output.split())
                self.assertEqual((ns, black, white, pairs, vertical),
                                 (lx*ly*lz, lx*ly*lz//2, lx*ly*lz//2,
                                  lx*ly*lz, lx*ly*lz))
        for lx, ly, lz in [(3, 4, 2), (4, 5, 2), (1, 4, 2), (4, 4, 0)]:
            with self.subTest(invalid_lattice=(lx, ly, lz)):
                result = subprocess.run(
                    [str(self.exe), "geometry", str(lx), str(ly), str(lz)],
                    text=True, capture_output=True
                )
                self.assertNotEqual(result.returncode, 0)

    def test_exchange_detailed_balance_from_full_padded_sse_weights(self):
        for beta_a, beta_b, na, nb in [(.5, 2., 3, 7), (2., .5, 7, 3), (.3, 3., 5, 5)]:
            m, shift = 12, 13.7
            def log_weight(beta, n, product):
                return -beta*shift+n*math.log(beta)+math.lgamma(m-n+1)-math.lgamma(m+1)+math.log(product)
            before = log_weight(beta_a, na, 2.3)+log_weight(beta_b, nb, 4.1)
            after = log_weight(beta_a, nb, 4.1)+log_weight(beta_b, na, 2.3)
            actual = float(subprocess.check_output([str(self.exe), 'swap', str(beta_a), str(beta_b), str(na), str(nb)], text=True))
            self.assertAlmostEqual(actual, after-before, places=12)
            self.assertAlmostEqual(before+min(0., actual), after+min(0., -actual), places=12)


if __name__ == "__main__": unittest.main()
