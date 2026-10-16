# Repository audit

The folder inventory enumerates all project folders, including historical sweep
subdirectories. Git internals, dependency environments, generated frontend distribution,
and rebuildable caches are excluded from per-folder enumeration.
See `folder-inventory.json` and `cleanup-manifest.json` for the inventory and exact
file moves/removals. Removed duplicate manuscript files were compared byte for byte;
removed build artifacts have SHA-256 hashes in the manifest.

## Decisions

- Preserve active C++/API/frontend code and all pre-existing uncommitted changes.
- Consolidate the two manuscript trees into `docs/manuscript`.
- Preserve prior kernels under `archive/legacy-kernels`; the active build does not use them.
- Preserve sweep inputs and historical scientific results under `backend/files`.
  These are not proven stale merely because they are old. Existing notebooks and
  the historical runner retain their original paths and conventions.
- Preserve API caches and temporary-run results. They are local output, not trusted
  benchmark references. New `.gitignore` rules prevent newly generated artifacts
  from being added; previously tracked generated artifacts, runtime caches and local configuration
  are removed from the index while available local copies are preserved.
- Preserve changed logs and output snapshots in ignored `data/legacy` rather than
  deleting possible research evidence.
- Retain Windows runtime DLLs beside the active binaries and retain editor settings.
- Move developer inputs and scripts into named folders; keep the old comparison
  command as a compatibility wrapper.
- Fix the visualization's hard-coded personal path, anchor sweep generation to its
  script directory, fix generated source paths and scheduler executable/arguments,
  and add active-build header dependencies and platform-specific linker flags.

## Benchmark readiness

Repository organization is separate from scientific validation. The benchmark guide
records required convergence, seed, uncertainty and cache controls. The existing
comparison script is preserved as a diagnostic report. No claim of SLQMC/ED agreement
or performance superiority follows from this cleanup.

## Validation performed

- Normal and SSE_DEBUG C++ binaries rebuilt successfully with GNU Make/g++.
- Direct 4x4 diagnostic run completed 200 equilibration and 200 measurement steps
  with exit code zero; output is preserved locally in `data/validation`.
- Project Python sources passed syntax parsing; all recorded move destinations exist.
- API import could not complete because the available Python lacks pandas.
- Frontend build could not run because npm dependencies are absent (`tsc` unavailable).
  Install the locked frontend dependencies with `npm ci` before rebuilding.
- Full API, UI and scientific agreement validation remain pending dependency setup.
