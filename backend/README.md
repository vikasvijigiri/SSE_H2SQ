# Backend

The C++ SSE kernel lives in `src/` and `include/`; the FastAPI service lives in `app/`.
ED reference code and research notebooks live in `Analysis/`.
See [the repository guide](../README.md) for build and run commands and
[the benchmark guide](../benchmarks/README.md) before collecting comparison results.

`files/beta_*` contains preserved historical sweep configurations and results.
The script `build/launch_dir.py` regenerates sweep inputs; running it overwrites matching inputs.
The older copied Makefiles and notebook runner use historical executable conventions.
Use the canonical build for new experiments.
