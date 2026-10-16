# QMC versus exact diagonalization

Run the API on port 8001, then `python benchmarks/compare_qmc_ed.py`.
The historical `python run_tests.py` command remains supported.
Install `requests` separately for this HTTP client: `python -m pip install requests`.

The existing script is a manual comparison report, not a statistical acceptance test.
It compares three 4x4 single-layer parameter sets with J3=0. The code calls the
simulation QMC/SSE; confirm which SLQMC algorithm/configuration is intended before
labeling published results.

## Protocol for the forthcoming benchmark

1. Rebuild the simulation and record the Git revision, compiler/flags, Python package
   versions, hardware, operating system, and complete input parameters.
2. Match lattice geometry and boundaries, beta, coupling conventions, field (hx/h),
   layer count (Lz/Nl), observable definitions and normalization. ED has no J3 term;
   restrict matched cases to J3=0. The API limits ED to 16 spins.
3. Establish the ED reference accuracy. Sparse calculations retain a limited set of
   eigenpairs; convergence in retained states is required at finite temperature.
   Do not treat a truncated thermal trace as automatically exact.
4. Use fresh computations rather than historical API caches. Preserve the cache files
   before clearing them. Record stochastic seeds and use independent runs; first
   confirm the current implementation exposes sufficient seed control.
5. Separate equilibration from measurement. Estimate autocorrelation and effective
   sample size; report uncertainty using blocking or another justified estimator.
   Check observable differences against combined uncertainties and predefined tolerances.
6. Measure compilation, ED construction/diagonalization, sampling, and end-to-end time
   separately. Report repeated-run statistics and approximation settings.
7. Save immutable raw output, input files and metadata in a unique directory under
   `benchmarks/results/`; promote selected validated summaries to versioned artifacts.

Percentage agreement in the current report is descriptive and unstable near zero.
The script does not establish convergence, uncertainty, reproducibility, or a pass/fail
criterion. No new scientific benchmark was performed during repository cleanup.
