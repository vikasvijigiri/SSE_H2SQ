# Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
"""Persistent, cancellable ED/QMC sweeps; numerical artifacts stay local."""
import asyncio
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import threading
import time
import uuid
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import scipy
from pydantic import BaseModel, Field, model_validator

from .physics_model import exact_reference, MODEL_VERSION
from .benchmark_statistics import summarize

ROOT = Path(__file__).resolve().parents[3]
STORE = ROOT / ".cache" / "benchmarks"
EXECUTABLE = ROOT / "backend/build" / ("sse_benchmark.exe" if os.name == "nt" else "sse_benchmark")
ACTIVE = {}
LOCK = threading.RLock()
METHODS = {"vertex": "Per-vertex pairing", "premarked": "Spatial premarking", "canonical": "Premarking + worldline Metropolis", "tempered": "Tempered premarking + worldline Metropolis"}
OBSERVABLES = {
    "E": ("Energy per spin", "E/N"), "E2": ("Second energy moment", "<(H/N)^2>"), "E4": ("Fourth energy moment", "<(H/N)^4>"),
    "mz2": ("Staggered magnetization squared", "<m_s^2>"), "mz4": ("Fourth staggered moment", "<m_s^4>"), "mx": ("Transverse magnetization", "<sigma_x>"),
    "cv": ("Heat capacity per spin", "C_V/N"), "stripe2": ("Stripe order squared", "<m_(pi,0)^2+m_(0,pi)^2>"),
    "polar2": ("Molecular polarization squared", "<|P|^2>"), "binder": ("Stripe Binder parameter", "Q"),
    "plaquette": ("Minimal plaquette loop", "<A_p>"), "ice": ("Ice-rule fraction", "f_ice"), "rho": ("Normalized ice-rule parameter", "rho/N_p"),
    "chi_stag": ("Staggered Kubo susceptibility", "chi_s/N"), "chi_polar": ("Dipole Kubo susceptibility", "chi_P/N_p"), "chi_x": ("Transverse susceptibility", "d<m_x>/dh"),
    **{f"ct{i}": (f"Dipole correlation at tau/beta={fraction:g}", "c(tau)") for i, fraction in enumerate([0, .125, .25, .375, .5])},
}


class BenchmarkRequest(BaseModel):
    Lx: int = 4
    Ly: int = 2
    J0: float = 1
    beta: float = Field(default=2, gt=0, le=100)
    J1: float = .3
    J2: float = .025
    hx: float = .2
    axes: dict[str, list[float]] = Field(default_factory=lambda: {"beta": [.5, 2, 8], "J1": [0, .15, .3], "J2": [0, .025, .05], "hx": [0, .2, .5]})
    algorithms: list[str] = Field(default_factory=lambda: ["vertex", "premarked", "canonical", "tempered"])
    chains: int = Field(default=4, ge=2, le=8)
    n_therm: int = Field(default=50000, ge=100, le=1000000)
    n_measure: int = Field(default=105000, ge=1000, le=1000000)
    seed: int = Field(default=20261016, ge=0, le=2147480000)

    @model_validator(mode="after")
    def validate_grid(self):
        if self.Lx < 2 or self.Ly < 2 or self.Lx % 2 or self.Ly % 2 or self.Lx*self.Ly > 16:
            raise ValueError("Use an even periodic lattice with at most 16 spins.")
        if not self.algorithms or any(method not in METHODS for method in self.algorithms) or len(set(self.algorithms)) != len(self.algorithms):
            raise ValueError("Choose distinct supported update schedules.")
        if not self.axes or any(axis not in {"beta", "J1", "J2", "hx"} for axis in self.axes):
            raise ValueError("Sweep beta, J1, J2 or hx.")
        if any(not values or len(values) > 12 or any(not np.isfinite(x) or abs(x) > 100 for x in values) for values in self.axes.values()) or sum(map(len, self.axes.values())) > 40:
            raise ValueError("Provide 1-12 finite values per axis, at most 40 points total.")
        if any(x <= 0 for x in self.axes.get("beta", [])):
            raise ValueError("Beta must be positive.")
        if self.Lx*self.Ly > 10 and (self.hx != 0 or any(x != 0 for x in self.axes.get("hx", []))):
            raise ValueError("Complete quantum ED is limited to 10 spins. Use 4x2 for field sweeps, or h=0 and omit the hx sweep on 4x4.")
        if not all(np.isfinite(x) for x in (self.J0, self.J1, self.J2, self.hx)):
            raise ValueError("Couplings must be finite.")
        return self


def save_job(job):
    with LOCK:
        directory = STORE / job["id"]
        directory.mkdir(parents=True, exist_ok=True)
        temporary = directory / "job.tmp"
        temporary.write_text(json.dumps(job, allow_nan=False), encoding="utf-8")
        # Windows readers/virus scanners may briefly deny a destination rename.
        for attempt in range(20):
            try:
                temporary.replace(directory / "job.json")
                break
            except PermissionError:
                if attempt == 19:
                    raise
                time.sleep(.05)


def get_job(identity):
    if not identity or any(c not in "0123456789abcdef" for c in identity) or len(identity) != 32:
        raise FileNotFoundError(identity)
    with LOCK:
        job = json.loads((STORE / identity / "job.json").read_text(encoding="utf-8"))
        if job["status"] in {"running", "queued"} and identity not in ACTIVE:
            job["status"] = "interrupted"
            job["message"] = "The service restarted; completed point results remain available."
            save_job(job)
    return job


def list_jobs():
    rows = []
    for path in STORE.glob("*/job.json"):
        job = get_job(path.parent.name)
        rows.append({key: job[key] for key in ("id", "created", "status", "request", "completed_runs", "total_runs")})
    return sorted(rows, key=lambda row: row["created"], reverse=True)[:20]


def provenance():
    sources = [ROOT / "backend/src/updates.cpp", ROOT / "backend/src/variables.cpp", ROOT / "backend/src/observables.cpp", ROOT / "backend/benchmarks/runner.cpp", Path(__file__).with_name("physics_model.py"), Path(__file__).with_name("benchmark_statistics.py")]
    sources += sorted((ROOT / "backend/include").glob("*.hpp"))
    return {"model_version": MODEL_VERSION, "revision": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(), "working_tree_dirty": bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT)), "engine_sha256": hashlib.sha256(EXECUTABLE.read_bytes()).hexdigest(), "source_sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources}, "python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__, "platform": platform.platform(), "processor": platform.processor(), "compiler": subprocess.check_output(["g++", "--version"], text=True).splitlines()[0], "flags": "-O2 -std=c++17 -Wall; no fast-math"}


def archive_sources(job):
    folder = STORE / job["id"]
    for name, digest in job["provenance"]["source_sha256"].items():
        data = (ROOT / name).read_bytes()
        if hashlib.sha256(data).hexdigest() != digest:
            raise RuntimeError("Source changed during experiment setup; deploy and retry.")
        path = folder / "sources" / digest / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
    shutil.copyfile(EXECUTABLE, folder / ("engine.exe" if os.name == "nt" else "engine"))


def run_chain(params, algorithm, seed, directory, request, cancelled):
    directory = directory.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    entries = {**params, "Lz": 1, "J3": 0, "n_therm": request.n_therm, "n_measure": request.n_measure}
    input_text = "\n".join(f"{key} {value}" for key, value in entries.items())+"\n"
    if all((directory / name).exists() for name in ("input.in", "metadata.json", "trace.npz")):
        metadata = json.loads((directory / "metadata.json").read_text())
        if metadata.get("seed") == seed and (directory / "input.in").read_text() == input_text:
            with np.load(directory / "trace.npz", allow_pickle=False) as data:
                trace = data["trace"]
            if len(trace) == request.n_measure and all(np.isfinite(trace[key]).all() for key in trace.dtype.names):
                return trace, metadata
    (directory / "input.in").write_text(input_text)
    started = time.monotonic()
    captured_engine = directory.parents[2] / ("engine.exe" if os.name == "nt" else "engine")
    process = subprocess.Popen([str(captured_engine if captured_engine.exists() else EXECUTABLE), "input.in", "trace.csv", algorithm, str(seed)], cwd=directory, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        while True:
            try:
                stdout, stderr = process.communicate(timeout=1)
                break
            except subprocess.TimeoutExpired:
                if cancelled.is_set():
                    process.kill()
                    process.communicate()
                    raise InterruptedError("Benchmark cancelled.")
                if time.monotonic()-started > 1800:
                    process.kill()
                    process.communicate()
                    raise RuntimeError("A QMC chain exceeded its 30-minute time limit.")
        if process.returncode:
            raise RuntimeError(stderr or stdout or "QMC process failed.")
    finally:
        if process.poll() is None:
            process.kill()
            process.wait()
    metadata = json.loads(stdout)
    metadata.update(seed=seed, wall_seconds=time.monotonic()-started)
    (directory / "metadata.json").write_text(json.dumps(metadata, indent=2))
    trace = pd.read_csv(directory / "trace.csv").to_records(index=False)
    if len(trace) != request.n_measure or not all(np.isfinite(trace[key]).all() for key in trace.dtype.names):
        raise RuntimeError("Incomplete or nonfinite QMC measurement trace.")
    np.savez_compressed(directory / "trace.npz", trace=trace)
    (directory / "trace.csv").unlink()
    return trace, metadata


def execute(job, request, cancelled):
    try:
        job["status"] = "running"
        save_job(job)
        base = {key: getattr(request, key) for key in ("Lx", "Ly", "J0", "beta", "J1", "J2", "hx")}
        reused = {}
        for axis, axis_values in request.axes.items():
            for value in sorted(set(axis_values)):
                if cancelled.is_set():
                    raise InterruptedError("Benchmark cancelled.")
                params = {**base, axis: value}
                identity = json.dumps({k: float(v) for k, v in params.items()}, sort_keys=True)
                if identity in reused:
                    job["points"].append({**reused[identity], "axis": axis, "x": value, "reused": True})
                    save_job(job)
                    continue
                point_id = hashlib.sha256(identity.encode()).hexdigest()[:16]
                point = {"axis": axis, "x": value, "params": params, "point_id": point_id, "ed": exact_reference(params), "methods": {}}
                job["points"].append(point)
                for algorithm in request.algorithms:
                    chains, metadata = [], []
                    job["message"] = f"{axis}={value:g} | {METHODS[algorithm]} | {request.chains} independent chains"
                    save_job(job)
                    completed = {}
                    with ThreadPoolExecutor(max_workers=min(4, request.chains)) as pool:
                        futures = {pool.submit(run_chain, params, algorithm,
                                   request.seed+chain*1009+int(point_id[:6], 16),
                                   STORE / job["id"] / point_id / algorithm / str(chain),
                                   request, cancelled): chain for chain in range(request.chains)}
                        try:
                            for future in as_completed(futures):
                                completed[futures[future]] = future.result()
                                job["completed_runs"] += 1
                                save_job(job)
                        except Exception:
                            cancelled.set()
                            raise
                    chains, metadata = zip(*(completed[chain] for chain in range(request.chains)))
                    summary = summarize(chains, params)
                    if any(row["cutoff_saturations"] for row in metadata):
                        for observation in summary["observables"].values():
                            observation.setdefault("flags", []).append("operator cutoff reached: rerun with more thermalization/headroom")
                            observation["status"] = "review"
                    if algorithm == "tempered" and any(row.get("replicas", 1) > 1 and row.get("round_trips", 0) < 10 for row in metadata):
                        for observation in summary["observables"].values():
                            observation.setdefault("flags", []).append("fewer than ten measured temperature round trips in a chain")
                            observation["status"] = "review"
                    summary["runs"] = metadata
                    summary["measurement_seconds"] = sum(row["measurement_seconds"] for row in metadata)
                    summary["total_seconds"] = sum(row["wall_seconds"] for row in metadata)
                    summary["update_seconds"] = {stage: sum(row["update_seconds"][stage] for row in metadata) for stage in ("diagonal", "cluster", "gauge", "canonical")}
                    compare_reference(summary, point["ed"])
                    point["methods"][algorithm] = summary
                    save_job(job)
                reused[identity] = point
        job.update(status="complete", message="All requested complete ED references and QMC chains finished.")
    except InterruptedError as error:
        job.update(status="cancelled", message=str(error))
    except Exception as error:
        job.update(status="failed", message=str(error))
    finally:
        job["finished"] = datetime.now(timezone.utc).isoformat()
        try:
            save_job(job)
        finally:
            with LOCK:
                ACTIVE.pop(job["id"], None)


async def start_job(request):
    if not EXECUTABLE.exists():
        raise RuntimeError("Build the benchmark engine with DEPLOY ALL first.")
    kernel_sources = [ROOT / "backend/benchmarks/runner.cpp", *(ROOT / "backend/src").glob("*.cpp"), *(ROOT / "backend/include").glob("*.hpp")]
    if any(path.stat().st_mtime > EXECUTABLE.stat().st_mtime for path in kernel_sources):
        raise RuntimeError("C++ source changed. Run DEPLOY ALL before starting an experiment.")
    run_provenance = provenance()
    with LOCK:
        if ACTIVE:
            raise RuntimeError("A benchmark is already running. Finish or cancel it before starting another.")
        if not EXECUTABLE.exists():
            raise RuntimeError("Build the benchmark engine with make -C backend/build benchmark.")
        identity = uuid.uuid4().hex
        event = threading.Event()
        ACTIVE[identity] = event
    unique = {json.dumps({key: float(value if key == axis else getattr(request, key)) for key in ("Lx", "Ly", "J0", "beta", "J1", "J2", "hx")}, sort_keys=True) for axis, values in request.axes.items() for value in values}
    job = {"id": identity, "created": datetime.now(timezone.utc).isoformat(), "status": "queued", "request": request.model_dump(), "completed_runs": 0, "total_runs": len(unique)*len(request.algorithms)*request.chains, "points": [], "message": "Queued", "provenance": run_provenance, "methods": METHODS, "observables": OBSERVABLES, "citations": [{"title": "Biswas & Damle: plaquette clusters and premarking", "url": "https://arxiv.org/abs/1812.05326"}, {"title": "Hukushima & Nemoto: exchange Monte Carlo", "url": "https://arxiv.org/abs/cond-mat/9512035"}, {"title": "Wolff: correlated-sample error analysis", "url": "https://arxiv.org/abs/hep-lat/0306017"}], "definitions": {"hamiltonian": "-J0 sum A + J1 sum(s0*s2+s1*s3) + J2 sum tilde(Pa).tilde(Pb) - hx sum sigma_x", "dipoles": "Cyclic (s0,s1,s2,s3); tilde Px=(s0-s1-s2+s3)/4, tilde Py=(s0+s1-s2-s3)/4; Pp=(-1)^x tilde Pp. Neighbor signs are opposite, so -J2 Pa.Pb equals +J2 tilde Pa.tilde Pb.", "ed_errors": "Complete finite-system ED has zero sampling error. Eigenpair residual is reported separately; no invented stochastic ED error bars.", "qmc_errors": "One standard error; maximum of correlation-aware, block-jackknife and independent-chain estimates. Review flagged points before publication.", "tempering": "Geometric temperature ladder from beta<=0.25/max(1,|couplings|) to target, adjacent ratio<=1.6. Equal padding lengths frozen after thermalization; exchange acceptance=min(1,(beta_i/beta_j)^(n_j-n_i)). Odd/even exchanges every ten sweeps. Measured round trips and per-edge acceptance are recorded; tune the ladder if overlap is poor.", "autocorrelation": "Complete update schedules are stationary kernels. Individual diagonal/cluster steps are timed, not assigned misleading isolated-chain autocorrelation times.", "stripe": "m_q=sum sigma_z exp(i q.r)/(2N); stripe2=m_(pi,0)^2+m_(0,pi)^2; Q=0.5[1-<stripe2^2>/(3<stripe2>^2)].", "ice": "rho/Np=(4 f_ice-1)/3", "susceptibility": "Diagonal Kubo responses use imaginary-time averages; equal-time beta*variance is not substituted at nonzero field.", "scope": "One layer, J3=0, periodic even checkerboard. One-factor-at-a-time sweeps; small-system validation does not establish thermodynamic phase boundaries."}}
    try:
        archive_sources(job)
        save_job(job)
    except Exception:
        with LOCK:
            ACTIVE.pop(identity, None)
        raise
    asyncio.create_task(asyncio.to_thread(execute, job, request, event))
    return job


async def resume_job(identity):
    job = get_job(identity)
    current = provenance()
    # An analysis-only change permits reanalysis; changing the physical model or
    # sampler requires a fresh experiment. Record both analysis versions.
    physical_sources = [key for key in current["source_sha256"] if not key.endswith("benchmark_statistics.py")]
    if any(job["provenance"]["source_sha256"].get(key) != current["source_sha256"][key] for key in physical_sources):
        raise RuntimeError("Sampler/reference source changed. Start a fresh experiment.")
    request = BenchmarkRequest(**job["request"])
    with LOCK:
        if ACTIVE:
            raise RuntimeError("Another benchmark is active.")
        event = threading.Event()
        ACTIVE[identity] = event
    job.setdefault("previous_provenance", []).append(job["provenance"])
    job.update(status="queued", message="Resuming verified stored chains and recomputing diagnostics", points=[], completed_runs=0, provenance=current)
    job.pop("finished", None)
    try:
        archive_sources(job)
        save_job(job)
    except Exception:
        with LOCK:
            ACTIVE.pop(identity, None)
        raise
    asyncio.create_task(asyncio.to_thread(execute, job, request, event))
    return job


def cancel_job(identity):
    with LOCK:
        event = ACTIVE.get(identity)
        if event:
            event.set()
    return bool(event)


def compare_reference(summary, reference):
    """Shared-sector traps can pass R-hat; an available exact reference adds a check."""
    for key, observation in summary["observables"].items():
        target = reference["values"].get(key)
        error = observation.get("error")
        observation["difference"] = observation["mean"]-target if observation["mean"] is not None and target is not None else None
        observation["pull"] = observation["difference"]/error if error is not None and error > 0 and observation["difference"] is not None else None
        observation["ess_per_second"] = observation.get("ess")/summary["total_seconds"] if observation.get("ess") is not None else None
        if observation["pull"] is not None and abs(observation["pull"]) > 4:
            flag = "ED discrepancy exceeds 4 SE: sampling/systematic review required"
            if flag not in observation.setdefault("flags", []):
                observation["flags"].append(flag)
            observation["status"] = "review"
    summary["validation_status"] = "review" if any(value["status"] == "review" for value in summary["observables"].values()) else "checks passed"


def csv_result(job):
    output = io.StringIO()
    writer = csv.writer(output)
    writer.writerow(["axis", "x", "Lx", "Ly", "beta", "J0", "J1", "J2", "hx", "algorithm", "observable", "ED", "ED_sampling_SE", "QMC", "QMC_SE", "tau_int_sweeps", "ESS", "Rhat", "status", "tau_SE_estimate", "ESS_per_summed_chain_second", "flags"])
    for point in job["points"]:
        for method, summary in point["methods"].items():
            for key, value in summary["observables"].items():
                writer.writerow([point["axis"], point["x"], *[point["params"][name] for name in ("Lx", "Ly", "beta", "J0", "J1", "J2", "hx")], method, key, point["ed"]["values"].get(key), 0, value["mean"], value["error"], value.get("tau_int"), value.get("ess"), value.get("rhat"), value["status"], value.get("tau_error"), value.get("ess_per_second"), "; ".join(value.get("flags", []))])
    return output.getvalue()


def figure_result(job, axis, observable, format, quantity="mean"):
    """Standalone vector figures with explicit physical normalization and 1-SE bars."""
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    import matplotlib as mpl
    mpl.rcParams.update({"pdf.fonttype": 42, "svg.fonttype": "none", "font.family": "serif", "mathtext.fontset": "stix"})
    points = sorted((p for p in job["points"] if p["axis"] == axis), key=lambda p: p["x"])
    figure = Figure(figsize=(6.4, 4.4), layout="constrained")
    figure.get_layout_engine().set(rect=(0, .065, 1, .935))
    FigureCanvasAgg(figure)
    axes = figure.subplots()
    axes.spines[["top", "right"]].set_visible(False)
    axes.grid(alpha=.15)
    if quantity == "mean":
        axes.plot([p["x"] for p in points], [p["ed"]["values"].get(observable) for p in points], "kD--", markersize=4, label="Complete ED (sampling SE = 0)")
    elif quantity == "difference":
        axes.axhline(0, color="black", linestyle="--", linewidth=1)
    colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]
    for method, color, marker in zip(METHODS, colors, ("o", "s", "^", "v")):
        rows = [(p["x"], p["methods"][method]["observables"].get(observable)) for p in points if method in p["methods"]]
        for x, row in rows:
            if row and row.get(quantity) is not None:
                error = row.get("error") if quantity in {"mean", "difference"} else row.get("tau_error") if quantity == "tau_int" else None
                axes.errorbar(x, row[quantity], yerr=error, fmt=marker if row["status"] != "review" else "x", color=color, capsize=3, markersize=5)
        if rows:
            axes.plot([], [], marker, color=color, label=METHODS[method])
    axes.set_xlabel({"beta": r"$\beta J_0$", "J1": r"$J_1$", "J2": r"$J_2$", "hx": r"$h_x$"}[axis] if axis != "beta" or job["request"]["J0"] == 1 else r"$\beta$")
    axes.set_ylabel(OBSERVABLES[observable][1] if quantity == "mean" else "QMC - ED" if quantity == "difference" else "Integrated correlation time (sweeps)" if quantity == "tau_int" else "ESS / summed chain second")
    fixed = {k: job["request"][k] for k in ("beta", "J0", "J1", "J2", "hx") if k != axis}
    axes.set_title(f'{job["request"]["Lx"]} x {job["request"]["Ly"]} periodic | '+", ".join(f"{k}={v:g}" for k, v in fixed.items()), fontsize=9)
    axes.legend(frameon=False, fontsize=7)
    figure.text(.02, .015, "Bars: 1 SE (tau: asymptotic estimate). Crosses require review; missing bars mean uncertainty unresolved.", fontsize=6)
    output = io.BytesIO()
    figure.savefig(output, format=format, dpi=300, metadata={"Creator": "SSE H2SQ benchmark"} if format != "png" else None)
    return output.getvalue()
