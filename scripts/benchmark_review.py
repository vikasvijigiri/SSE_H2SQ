# Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
"""Apply current exact-reference diagnostics to a completed local experiment."""
import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "backend"))
from app.services.benchmark import get_job, save_job, compare_reference


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("identity")
    args = parser.parse_args()
    job = get_job(args.identity)
    if job["status"] != "complete":
        raise SystemExit("Finish the experiment before applying an offline review.")
    checked, flagged = 0, 0
    for point in job["points"]:
        for method in point["methods"].values():
            compare_reference(method, point["ed"])
            for value in method["observables"].values():
                checked += 1
                flagged += value["status"] == "review"
    job["review"] = {"ed_discrepancy_threshold_SE": 4, "checked": checked, "flagged": flagged,
                     "note": "Reapplied exact-reference comparisons to stored means and uncertainties; no numerical samples changed."}
    save_job(job)
    print(f"Reviewed {checked} observable comparisons; {flagged} require review.")


if __name__ == "__main__": main()
