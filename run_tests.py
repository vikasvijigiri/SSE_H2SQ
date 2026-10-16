"""Compatibility entry point for the QMC/ED comparison script."""
from pathlib import Path
import runpy

if __name__ == "__main__":
    runpy.run_path(str(Path(__file__).resolve().parent / "benchmarks" / "compare_qmc_ed.py"), run_name="__main__")
