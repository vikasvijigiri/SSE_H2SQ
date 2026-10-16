"""Run with python -m unittest discover -s scripts/tests."""
import math
from pathlib import Path
import sys
import tempfile
import unittest
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))
from app.services import qmc_results


class QmcResultsTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.previous = qmc_results.DATABASE
        qmc_results.DATABASE = Path(self.temporary.name) / "runs.sqlite3"

    def tearDown(self):
        qmc_results.DATABASE = self.previous
        self.temporary.cleanup()

    def test_jackknife_matches_explicit_leave_one_out(self):
        samples = [1.0, 3.0, 5.0]
        for value in samples:
            qmc_results.save_run("same-physics", {"enrg": value})
        qmc_results.save_run("other-physics", {"enrg": 1000})
        estimate = qmc_results.jackknife_summary("same-physics")
        leave_one_out = [(sum(samples) - value) / 2 for value in samples]
        center = sum(leave_one_out) / 3
        expected_error = math.sqrt(2 / 3 * sum((value - center) ** 2 for value in leave_one_out))
        self.assertEqual(estimate["n_runs"], 3)
        self.assertEqual(estimate["observables"]["enrg"]["mean"], 3)
        self.assertAlmostEqual(estimate["observables"]["enrg"]["standard_error"], expected_error)

    def test_legacy_import_deduplicates_and_single_run_has_no_error(self):
        record = {"enrg": -0.8, "timestamp": "2026-10-16T12:00:00"}
        qmc_results.import_legacy("key", record)
        qmc_results.import_legacy("key", record)
        estimate = qmc_results.jackknife_summary("key")
        self.assertEqual(estimate["n_runs"], 1)
        self.assertIsNone(estimate["observables"]["enrg"]["standard_error"])

    def test_raw_file_and_concurrent_runs_are_preserved(self):
        def store(value):
            return qmc_results.save_run("key", {"enrg": value}, f"E\n{value}\n")
        with ThreadPoolExecutor(max_workers=4) as pool:
            records = list(pool.map(store, range(8)))
        self.assertEqual(qmc_results.jackknife_summary("key")["n_runs"], 8)
        for record in records:
            self.assertEqual(qmc_results.result_text(record["run_id"]), f"E\n{record['enrg']}\n")
            qmc_results.import_legacy("key", record)
        self.assertEqual(qmc_results.jackknife_summary("key")["n_runs"], 8)
        self.assertIsNone(qmc_results.result_text("../../unrelated-file"))


if __name__ == "__main__":
    unittest.main()
