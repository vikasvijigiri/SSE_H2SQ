# SSE-H2SQ

SSE Quantum Monte Carlo simulation, exact diagonalization reference calculations,
and a web interface for the H2SQ model.

## Repository layout

| Directory | Purpose |
| --- | --- |
| `backend/src`, `backend/include` | Active C++ simulation kernel |
| `backend/app` | FastAPI service and simulation orchestration |
| `backend/Analysis` | ED implementation and research notebooks; paths retained for notebook compatibility |
| `backend/build` | Canonical Makefile, sweep generator, scheduler template, local binaries |
| `backend/files` | Preserved sweep inputs, historical output, API caches and temporary runs |
| `frontend` | React/TypeScript interface |
| `benchmarks` | QMC-versus-ED comparison workflow |
| `configs/diagnostics` | Preserved simulation input examples |
| `scripts/diagnostics` | Manual Python simulation diagnostics |
| `tools/visualization` | Lattice figure generator and SVG |
| `docs/manuscript` | Canonical manuscript, bibliography, figures and references |
| `archive` | Historical kernels and site-specific scheduler helper |
| `data/legacy` | Local preserved output and logs, ignored by Git |

## Build and run

Requires Python with the packages in `backend/requirements.txt`, a C++17 compiler,
GNU Make, and Node.js/npm for the interface.

```sh
python -m pip install -r backend/requirements.txt
make -C backend/build
python -m uvicorn app.main:app --app-dir backend --port 8001
```

In another terminal:

```sh
cd frontend
npm ci
npm run dev
```

For a direct simulation, create a dedicated output directory and execute
`backend/build/sse_sim_test` (`.exe` on Windows) using absolute executable and input
paths. Results are written to the current working directory.
`make -C backend/build debugbin` builds the diagnostic binary.

## Validation and research

Read [the benchmark protocol](benchmarks/README.md) before interpreting QMC/ED comparisons.
`python run_tests.py` remains an alias for `python benchmarks/compare_qmc_ed.py`.
Both require the API to be running and execute potentially lengthy simulations.

See [the repository audit](docs/repository-audit.md) for cleanup scope, preserved
historical material, and outstanding benchmark requirements. Existing local physics and UI edits are included alongside the repository cleanup.

## Paper folder

The organized paper sources are in [docs/manuscript](docs/manuscript/README.md).
Open `main.tex` as the main document. Sections, figures, bibliography, and background
references have dedicated folders. Run `python scripts/validate_paper.py` to check
source paths and citations without creating an archive.
