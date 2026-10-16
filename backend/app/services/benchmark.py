# Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
"""Persistent, cancellable ED/QMC sweeps; numerical artifacts stay local."""
import asyncio
import csv
import hashlib
import io
import json
import math
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
from scipy.stats import poisson
from pydantic import BaseModel, Field, model_validator

from .physics_model import exact_reference, MODEL_VERSION
from .benchmark_statistics import summarize

ROOT = Path(__file__).resolve().parents[3]
STORE = ROOT / ".cache" / "benchmarks"
EXECUTABLE = ROOT / "backend/build" / ("sse_benchmark.exe" if os.name == "nt" else "sse_benchmark")
ACTIVE = {}
LOCK = threading.RLock()
METHODS = {"vertex": "Per-vertex pairing", "premarked": "Spatial premarking", "canonical": "Premarking + worldline Metropolis", "tempered": "Tempered premarking + worldline Metropolis", "symmetry":"Premarking + worldline + exact site-sublattice symmetry"}
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
    algorithms: list[str] = Field(default_factory=lambda: ["canonical","symmetry"])
    chains: int = Field(default=8, ge=2, le=16)
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
        if self.Lx*self.Ly > 10 and (self.Lx,self.Ly)!=(4,4) and (self.hx != 0 or any(x != 0 for x in self.axes.get("hx", []))):
            raise ValueError("Quantum references above 10 spins currently require the 4x4 symmetry-resolved lattice.")
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
    sources = [ROOT / "backend/src/updates.cpp", ROOT / "backend/src/variables.cpp", ROOT / "backend/src/observables.cpp", ROOT / "backend/benchmarks/runner.cpp", Path(__file__).with_name("physics_model.py"), Path(__file__).with_name("symmetry_ed.py"), Path(__file__).with_name("benchmark_statistics.py"), Path(__file__)]
    sources += sorted((ROOT / "backend/include").glob("*.hpp"))
    return {"model_version": MODEL_VERSION, "revision": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(), "working_tree_dirty": bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT)), "engine_sha256": hashlib.sha256(EXECUTABLE.read_bytes()).hexdigest(), "source_sha256": {path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources}, "python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__, "platform": platform.platform(), "processor": platform.processor(), "compiler": subprocess.check_output(["g++", "--version"], text=True).splitlines()[0], "flags": "-O2 -std=c++17 -Wall; no fast-math"}


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
                    if algorithm == "tempered" and any(any(n >= 100 and a/n < .02 for n, a in zip(row.get("swap_attempts", []), row.get("swap_accepted", []))) for row in metadata):
                        for observation in summary["observables"].values():
                            observation.setdefault("flags", []).append("temperature ladder bottleneck: exchange acceptance below 2 percent")
                            observation["status"] = "review"
                    summary["runs"] = metadata
                    summary["measurement_seconds"] = sum(row["measurement_seconds"] for row in metadata)
                    summary["total_seconds"] = sum(row["wall_seconds"] for row in metadata)
                    summary["update_seconds"] = {stage: sum(row["update_seconds"].get(stage,0.) for row in metadata) for stage in ("diagonal", "cluster", "gauge", "canonical","symmetry")}
                    sampling_checks(summary,metadata,algorithm)
                    compare_reference(summary, point["ed"])
                    if len(chains)>=8:
                        summary["error_scaling"]=precision_history(chains,params,metadata,algorithm,point["ed"],summary)
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
    job = {"id": identity, "created": datetime.now(timezone.utc).isoformat(), "status": "queued", "request": request.model_dump(), "completed_runs": 0, "total_runs": len(unique)*len(request.algorithms)*request.chains, "points": [], "message": "Queued", "provenance": run_provenance, "methods": METHODS, "observables": OBSERVABLES, "citations": [{"title": "Biswas & Damle: plaquette clusters and premarking", "url": "https://arxiv.org/abs/1812.05326"}, {"title": "Hukushima & Nemoto: exchange Monte Carlo", "url": "https://arxiv.org/abs/cond-mat/9512035"}, {"title": "Wolff: correlated-sample error analysis", "url": "https://arxiv.org/abs/hep-lat/0306017"}], "definitions": {"hamiltonian": "-J0 sum A + J1 sum(s0*s2+s1*s3) + J2 sum tilde(Pa).tilde(Pb) - hx sum sigma_x", "dipoles": "Cyclic (s0,s1,s2,s3); tilde Px=(s0-s1-s2+s3)/4, tilde Py=(s0+s1-s2-s3)/4; Pp=(-1)^x tilde Pp. Neighbor signs are opposite, so -J2 Pa.Pb equals +J2 tilde Pa.tilde Pb.", "ed_errors": "Complete finite-system ED has zero sampling error. Eigenpair residual is reported separately; no invented stochastic ED error bars.", "qmc_errors": "One standard error; maximum of correlation-aware, block-jackknife and independent-chain estimates. Review flagged points before publication.", "tempering": "Temperature ladder starts at beta<=0.25/max(1,|couplings|). Adjacent sqrt(beta) spacing<=0.7/sqrt(G_bound), also ratio<=1.6, with at most 64 replicas. G_bound is a conservative expansion-order density bound. Actual overlap and round trips must still be checked. Equal padding lengths frozen after thermalization; exchange acceptance=min(1,(beta_i/beta_j)^(n_j-n_i)). Odd/even exchanges every ten sweeps. Measured round trips and per-edge acceptance are recorded; tune the ladder if overlap is poor.", "autocorrelation": "Complete update schedules are stationary kernels. Individual diagonal/cluster steps are timed, not assigned misleading isolated-chain autocorrelation times.", "stripe": "m_q=sum sigma_z exp(i q.r)/(2N); stripe2=m_(pi,0)^2+m_(0,pi)^2; Q=0.5[1-<stripe2^2>/(3<stripe2>^2)].", "ice": "rho/Np=(4 f_ice-1)/3", "susceptibility": "Diagonal Kubo responses use imaginary-time averages; equal-time beta*variance is not substituted at nonzero field.", "scope": "One layer, J3=0, periodic even checkerboard. One-factor-at-a-time sweeps; small-system validation does not establish thermodynamic phase boundaries."}}
    job["citations"].extend([{"title":"Sandvik: computational studies and improved estimators","url":"https://arxiv.org/abs/1101.3281"},{"title":"Sandvik: SSE estimators, equations 24 and 27","url":"https://arxiv.org/abs/1909.10591"}])
    job["definitions"].update({"improved_estimators":"All nonidentity propagation indices are averaged for diagonal moments. Exact SSE Kubo sums replace random time embeddings. Dipole c(tau) integrates binomial time separations and cyclic origins analytically. Boundary estimates are retained for measured SE comparisons.","site_symmetry":"Independent flips of all even/odd SITE-parity spins generate an exact Z2 x Z2 symmetry, distinct from molecular A/B. J0 has two sites of each parity; retained J1/J2/J3 bonds join equal parities; transverse operators commute with the flips. Each of four orbit elements is chosen with probability 1/4, preserving every SSE matrix element and giving an acceptance-one reversible update. Runtime geometry checks and exhaustive shifted-weight tests verify its applicability. It supplements premarked clusters and worldline Metropolis; it is not a proof of rapid mixing between every inequivalent sector.","ed_errors":"ED has zero sampling error. 4x4 quantum references retain every eigenstate in the partition function and static observables. Lehmann responses retain all intermediate states linked to retained thermal endpoints; an exact-probability thermal-tail bound and eigenpair residual are reported separately.","scope":"One 2D periodic even checkerboard layer, J3=0. Dense complete ED up to 10 spins; full symmetry-resolved ED on 4x4. A finite nonzero-coupling grid does not validate a continuum of parameters or thermodynamic phase boundaries."})
    job["definitions"].update({"ergodicity": "For the certified schedule, each diagonal slot and whole-worldline proposal has a 1% identity probability. Any valid finite-padding configuration can reach an hx-only string by removing diagonal operators while skipping other slots; cluster flips can then make every worldline constant, and removal reaches the empty string. Free spins and the reverse moves reach every positive-weight configuration. Positive identity probabilities give aperiodicity. This is a constructive finite-support ergodicity argument, not a uniform mixing-time bound.", "SSE_cutoff": "For K=C-H=sum of positive semidefinite shifted operators, ||K||<=G. The expansion tail is bounded by P(Poisson(beta*G)>M). The certified schedule initializes M>=ceil(beta*G+12*sqrt(beta*G+1)+32), allows thermal growth, and freezes padding during measurements. Tail-probability and energy-moment bias bounds are recorded separately from sampling SE."})
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
    if job.get("source_experiments"):
        raise RuntimeError("A comparison retains original experiments; start a new experiment to resample.")
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


def combine_jobs(identities):
    """Compare completed schedules while retaining each source's data/provenance."""
    if not 2 <= len(identities) <= 4 or len(set(identities)) != len(identities):
        raise RuntimeError("Choose two to four distinct completed experiments.")
    jobs = [get_job(identity) for identity in identities]
    if any(job["status"] != "complete" for job in jobs):
        raise RuntimeError("Only completed experiments can be combined.")
    if any(len({job["request"][key] for job in jobs}) != 1 for key in ("chains", "n_therm", "n_measure", "seed")):
        raise RuntimeError("Use matching chain counts, sweep budgets and seed conventions for a controlled comparison.")
    if len({job["provenance"]["model_version"] for job in jobs}) != 1 or len({next((value for name, value in job["provenance"]["source_sha256"].items() if name.replace(chr(92), "/") == "backend/app/services/physics_model.py"), None) for job in jobs}) != 1:
        raise RuntimeError("Reference Hamiltonian versions differ.")
    def key(point):
        return json.dumps({name: float(value) for name, value in point["params"].items()}, sort_keys=True)
    lookups = [{key(point): point for point in job["points"]} for job in jobs]
    if any(set(lookup) != set(lookups[0]) for lookup in lookups):
        raise RuntimeError("The physical parameter grids differ; compare matching grids.")
    algorithms = [algorithm for job in jobs for algorithm in job["request"]["algorithms"]]
    if len(algorithms) != len(set(algorithms)):
        raise RuntimeError("Schedules overlap. Choose experiments with distinct schedules.")
    result = json.loads(json.dumps(jobs[0]))
    result["id"] = uuid.uuid4().hex
    result["created"] = datetime.now(timezone.utc).isoformat()
    result["request"]["algorithms"] = algorithms
    result["completed_runs"] = sum(job["completed_runs"] for job in jobs)
    result["total_runs"] = sum(job["total_runs"] for job in jobs)
    result["message"] = "Comparison of completed experiments; original traces and source provenance retained."
    result["source_experiments"] = [{"id": identity, "provenance": job["provenance"], "request": job["request"]} for identity, job in zip(identities, jobs)]
    result["definitions"].update(jobs[-1]["definitions"])
    result["citations"] = list({source["url"]: source for job in jobs for source in job["citations"]}.values())
    result.pop("review", None)
    for point in result["points"]:
        merged = {}
        for identity, lookup in zip(identities, lookups):
            source = lookup[key(point)]
            for name, reference in point["ed"]["values"].items():
                other = source["ed"]["values"].get(name)
                if reference is not None and (other is None or not np.isclose(reference, other, atol=1e-12, rtol=1e-12)):
                    raise RuntimeError("Complete ED references disagree.")
            for method, summary in source["methods"].items():
                merged[method] = {**summary, "source_experiment": summary.get("source_experiment", identity), "source_point": summary.get("source_point", source["point_id"])}
        point["methods"] = merged
    save_job(result)
    return result


def compare_reference(summary, reference):
    """Shared-sector traps can pass R-hat; an available exact reference adds a check."""
    for key, observation in summary["observables"].items():
        target = reference["values"].get(key)
        error = observation.get("error")
        observation["difference"] = observation["mean"]-target if observation["mean"] is not None and target is not None else None
        observation["pull"] = observation["difference"]/error if error is not None and error > 0 and observation["difference"] is not None else None
        observation["ess_per_second"] = observation.get("ess")/summary["total_seconds"] if observation.get("ess") is not None else None
        reference_bound=reference.get("absolute_error_bounds",{}).get(key,0.)
        sse_bound=summary.get("sse_moment_bias_bounds",{}).get(key)
        bound=reference_bound+(sse_bound or 0.)
        observation["reference_absolute_error_bound"]=reference_bound
        observation["sse_absolute_bias_bound"]=sse_bound
        if observation["pull"] is not None and abs(observation["difference"]) > 4*error+bound:
            flag = "ED discrepancy exceeds 4 SE: sampling/systematic review required"
            if flag not in observation.setdefault("flags", []):
                observation["flags"].append(flag)
            observation["status"] = "review"
    summary["validation_status"] = "review" if any(value["status"] == "review" for value in summary["observables"].values()) else "checks passed"


def trace_location(identity,point,method,chain):
    summary=point["methods"][method]
    run=summary["runs"][chain]
    return (run.get("source_experiment",summary.get("source_experiment",identity)),
            run.get("source_point",summary.get("source_point",point["point_id"])),
            int(run.get("source_chain",chain)))


def sampling_checks(summary,runs,algorithm):
    flags=[]
    cutoff_bounds=[]
    moment_bounds={key:[] for key in ("E","E2","E4")}
    for run in runs:
        if "operator_cutoff" not in run or "norm_mean_bound" not in run:continue
        cutoff,mu=int(run["operator_cutoff"]),float(run["norm_mean_bound"])
        tail=float(poisson.sf(cutoff,mu))
        run["sse_tail_probability_bound"]=tail if math.isfinite(tail) else None
        if not math.isfinite(tail):
            flags.append("invalid SSE cutoff-tail bound");continue
        cutoff_bounds.append(tail)
        if tail>1e-20:
            flags.append("SSE cutoff-tail bound exceeds 1e-20")
        if tail<1 and all(key in run for key in ("energy_shift_total","inverse_temperature","spins")):
            beta,n=run["inverse_temperature"],run["spins"]
            shift=abs(run["energy_shift_total"])
            delta={k:mu**k*(float(poisson.sf(cutoff-k,mu))+tail)/(1-tail) for k in range(1,5)}
            run["sse_moment_bias_bounds"]={}
            for power,name in ((1,"E"),(2,"E2"),(4,"E4")):
                value=sum(math.comb(power,k)*shift**(power-k)*delta[k]/beta**k for k in range(1,power+1))/n**power
                moment_bounds[name].append(value);run["sse_moment_bias_bounds"][name]=value
    summary["sse_cutoff_tail_bound"]=max(cutoff_bounds) if len(cutoff_bounds)==len(runs) else None
    summary["sse_moment_bias_bounds"]={name:max(values) for name,values in moment_bounds.items() if len(values)==len(runs)}
    if any(row["cutoff_saturations"] for row in runs):
        flags.append("operator cutoff reached: rerun with more thermalization/headroom")
    if algorithm=="tempered":
        if any(row.get("replicas",1)>1 and row.get("round_trips",0)<10 for row in runs):
            flags.append("fewer than ten measured temperature round trips in a chain")
        if any(any(n>=100 and a/n<.02 for n,a in zip(row.get("swap_attempts",[]),row.get("swap_accepted",[]))) for row in runs):
            flags.append("temperature ladder bottleneck: exchange acceptance below 2 percent")
    if algorithm=="symmetry":
        if any(not 0<row.get("local_skip_probability",0)<1 for row in runs):
            flags.append("positive local laziness missing from the certified schedule")
        if any(len(row.get("site_symmetry_orbit_counts",[]))!=4 or min(row.get("site_symmetry_orbit_counts",[0]))==0 for row in runs):
            flags.append("site-symmetry orbit activity missing")
        if sum(sum(row.get("site_symmetry_orbit_counts",[])) for row in runs)!=summary["retained_samples"]:
            flags.append("site-symmetry move count differs from measurement rows")
    for observation in summary["observables"].values():
        for flag in flags:
            if flag not in observation.setdefault("flags",[]):observation["flags"].append(flag)
        if flags:observation["status"]="review"
    if flags:summary["validation_status"]="review"


def precision_history(traces,params,runs,method,reference,full):
    history=[]
    for count in sorted(set([4,8,len(traces)])):
        if count<4 or count>len(traces):continue
        current=full if count==len(traces) else summarize(traces[:count],params)
        if count!=len(traces):
            current["total_seconds"]=sum(r["wall_seconds"] for r in runs[:count])
            sampling_checks(current,runs[:count],method)
            compare_reference(current,reference)
        history.append({"chains":count,"samples":sum(len(t) for t in traces[:count]),
            "observables":{name:{k:o.get(k) for k in ("mean","error","tau_int","ess","rhat","status")} for name,o in current["observables"].items()}})
    for row in history:
        row["expected_se_ratio"]=float(np.sqrt(history[0]["samples"]/row["samples"]))
        for name,o in row["observables"].items():
            original=history[0]["observables"][name]["error"]
            o["se_ratio"]=o["error"]/original if original and o["error"] is not None else None
    return history


def pool_jobs(identities):
    """Reanalyse original independent traces, never average reported SE bars."""
    if not 2<=len(identities)<=12 or len(set(identities))!=len(identities):
        raise RuntimeError("Choose two to twelve distinct completed experiments.")
    jobs=[get_job(identity) for identity in identities]
    if any(job["status"]!="complete" for job in jobs):
        raise RuntimeError("Only completed experiments can be pooled.")
    if len({job["provenance"]["model_version"] for job in jobs})!=1:
        raise RuntimeError("Hamiltonian versions differ; do not pool these experiments.")
    def key(point):
        return json.dumps({name:float(value) for name,value in point["params"].items()},sort_keys=True)
    lookups=[{key(point):point for point in job["points"]} for job in jobs]
    if any(set(lookup)!=set(lookups[0]) for lookup in lookups):
        raise RuntimeError("Choose identical physical parameter grids for pooling.")
    result=json.loads(json.dumps(jobs[0]))
    result.update(id=uuid.uuid4().hex,created=datetime.now(timezone.utc).isoformat(),
                  message="Pooled independent raw traces; all rows enter means, errors account for correlations.",
                  aggregation="independent-trace pooling",provenance=provenance())
    archive_sources(result)
    result["source_experiments"]=[{"id":j["id"],"provenance":j["provenance"],"request":j["request"]} for j in jobs]
    result["request"]["algorithms"]=list(dict.fromkeys(a for j in jobs for a in j["request"]["algorithms"]))
    completed,max_chains=0,0
    analysed={}
    for point in result["points"]:
        identity=key(point)
        if identity in analysed:
            point["methods"]=analysed[identity]
            continue
        sources=[(job,lookup[identity]) for job,lookup in zip(jobs,lookups)]
        for job,other in sources:
            for name,target in point["ed"]["values"].items():
                value=other["ed"]["values"].get(name)
                if target is not None and (value is None or not np.isclose(target,value,rtol=1e-10,atol=1e-10)):
                    raise RuntimeError("Exact references disagree; review model conventions before pooling.")
        methods={}
        for method in result["request"]["algorithms"]:
            traces,runs,seen=[],[],set()
            for job,source in sources:
                summary=source["methods"].get(method)
                if not summary:continue
                for chain,run in enumerate(summary["runs"]):
                    location=trace_location(job["id"],source,method,chain)
                    # A comparison may point at a trace already selected elsewhere.
                    if location in seen:
                        raise RuntimeError("A raw trace was selected twice through overlapping experiments.")
                    if any(previous["seed"]==run["seed"] for previous in runs):
                        raise RuntimeError("Repeated seed within a physical point/schedule; independence is not established.")
                    seen.add(location)
                    with np.load(STORE/location[0]/location[1]/method/str(location[2])/"trace.npz",allow_pickle=False) as data:
                        traces.append(data["trace"])
                    runs.append({**run,"source_experiment":location[0],"source_point":location[1],"source_chain":location[2]})
            # Older traces may not contain the optional boundary-estimator audit
            # columns. The shared primary estimator schema is required.
            common=[name for name in traces[0].dtype.names if all(name in trace.dtype.names for trace in traces)]
            traces=[trace[common] for trace in traces]
            summary=summarize(traces,point["params"])
            summary.update(runs=runs,measurement_seconds=sum(r["measurement_seconds"] for r in runs),
                           total_seconds=sum(r["wall_seconds"] for r in runs),
                           update_seconds={stage:sum(r["update_seconds"].get(stage,0.) for r in runs) for stage in ("diagonal","cluster","gauge","canonical","symmetry")})
            sampling_checks(summary,runs,method)
            compare_reference(summary,point["ed"])
            summary["error_scaling"]=precision_history(traces,point["params"],runs,method,point["ed"],summary)
            methods[method]=summary
            max_chains=max(max_chains,len(traces));completed+=len(traces)
        point["methods"]=methods;analysed[identity]=methods
    result["request"]["chains"]=max_chains
    result["completed_runs"]=result["total_runs"]=completed
    result["definitions"]["pooling"]="Raw traces are pooled within each physical point and update schedule. Duplicate traces/seeds are rejected. Means use all rows with sample-count weights; blocked jackknife, ACF and independent-chain estimates determine SE. Error reduction is measured, not forced to be monotone. Different update schedules remain separate."
    result["finished"]=datetime.now(timezone.utc).isoformat()
    save_job(result)
    return result


def csv_result(job):
    output = io.StringIO()
    writer = csv.writer(output)
    writer.writerow(["axis", "x", "Lx", "Ly", "beta", "J0", "J1", "J2", "hx", "algorithm", "observable", "ED", "ED_sampling_SE", "QMC", "QMC_SE", "tau_int_sweeps", "ESS", "Rhat", "status", "tau_SE_estimate", "ESS_per_summed_chain_second", "flags","ED_Lehmann_tail_absolute_bound","QMC_SSE_moment_absolute_bias_bound"])
    for point in job["points"]:
        for method, summary in point["methods"].items():
            for key, value in summary["observables"].items():
                writer.writerow([point["axis"], point["x"], *[point["params"][name] for name in ("Lx", "Ly", "beta", "J0", "J1", "J2", "hx")], method, key, point["ed"]["values"].get(key), 0, value["mean"], value["error"], value.get("tau_int"), value.get("ess"), value.get("rhat"), value["status"], value.get("tau_error"), value.get("ess_per_second"), "; ".join(value.get("flags", [])),point["ed"].get("absolute_error_bounds",{}).get(key,0.),value.get("sse_absolute_bias_bound")])
    return output.getvalue()


def overlay_controls(primary_identity,control_identities):
    """Keep a complete primary grid and attach available frozen control points."""
    identities=[primary_identity,*control_identities]
    if not control_identities or len(set(identities))!=len(identities):
        raise RuntimeError("Choose a primary study and distinct control experiments.")
    jobs=[get_job(identity) for identity in identities]
    primary=jobs[0]
    if primary["status"]!="complete" or any(job["status"] in {"running","queued"} for job in jobs[1:]):
        raise RuntimeError("Complete the primary study and finish or cancel each control first.")
    if len({job["provenance"]["model_version"] for job in jobs})!=1:
        raise RuntimeError("Hamiltonian conventions differ.")
    def key(point):return json.dumps({k:float(v) for k,v in point["params"].items()},sort_keys=True)
    lookups=[{key(point):point for point in job["points"]} for job in jobs]
    result=json.loads(json.dumps(primary));result.update(id=uuid.uuid4().hex,created=datetime.now(timezone.utc).isoformat(),
        aggregation="complete primary grid with available control summaries",
        message="Complete primary validation grid; original control points and their convergence warnings are retained.")
    result["source_experiments"]=[{"id":job["id"],"status":job["status"],"completed_runs":job["completed_runs"],"total_runs":job["total_runs"],"provenance":job["provenance"],"request":job["request"]} for job in jobs]
    locations=set();methods=set();maximum=0
    for point in result["points"]:
        merged={}
        for job,lookup in zip(jobs,lookups):
            source=lookup.get(key(point))
            if not source:continue
            for name,target in point["ed"]["values"].items():
                value=source["ed"]["values"].get(name)
                if target is not None and (value is None or not np.isclose(target,value,rtol=1e-10,atol=1e-10)):
                    raise RuntimeError("ED references disagree at a shared control point.")
            for method,summary in source["methods"].items():
                if method in merged:raise RuntimeError("Overlapping schedules require raw-trace pooling, not an overlay.")
                if not summary.get("runs"):continue
                merged[method]={**summary,"source_experiment":summary.get("source_experiment",job["id"]),"source_point":summary.get("source_point",source["point_id"])}
                for chain in range(len(summary["runs"])):
                    origin=trace_location(job["id"],source,method,chain)
                    locations.add((*origin,method))
                methods.add(method);maximum=max(maximum,len(summary["runs"]))
        point["methods"]=merged
    result["completed_runs"]=result["total_runs"]=len(locations)
    result["request"]["algorithms"]=[method for method in METHODS if method in methods]
    result["request"]["chains"]=maximum
    result["definitions"]["controls"]="The complete primary grid is retained. A stopped control contributes only already completed point summaries; missing controls are not invented. All source experiment statuses, original traces and convergence warnings remain available. Schedules are compared separately."
    save_job(result)
    return result


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
        if observable=="E" and any(p["ed"].get("ground_energy_per_spin") is not None for p in points):
            axes.plot([p["x"] for p in points],[p["ed"].get("ground_energy_per_spin") for p in points],color="#6b7280",linestyle=":",label="ED ground energy")
    elif quantity == "difference":
        axes.axhline(0, color="black", linestyle="--", linewidth=1)
    colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7","#E69F00"]
    for method, color, marker in zip(METHODS, colors, ("o", "s", "^", "v","P")):
        rows = [(p["x"], p["methods"][method]["observables"].get(observable)) for p in points if method in p["methods"]]
        for x, row in rows:
            if row and row.get(quantity) is not None:
                error = row.get("error") if quantity in {"mean", "difference"} else row.get("tau_error") if quantity == "tau_int" else None
                axes.errorbar(x, row[quantity], yerr=error, fmt=marker if row["status"] != "review" else "x", color=color, capsize=3, markersize=5)
        if rows:
            axes.plot([], [], marker, color=color, label=METHODS[method])
    axes.set_xlabel({"beta": r"$\beta J_0$", "J1": r"$J_1$", "J2": r"$J_2$", "hx": r"$h_x$"}[axis] if axis != "beta" or job["request"]["J0"] == 1 else r"$\beta$")
    axes.set_ylabel(OBSERVABLES[observable][1] if quantity == "mean" else "QMC - ED" if quantity == "difference" else "Integrated correlation time (sweeps)" if quantity == "tau_int" else "ESS / summed chain second")
    if quantity=="tau_int":
        positive=[summary["observables"].get(observable,{}).get("tau_int") for point in points for summary in point["methods"].values()]
        positive=[value for value in positive if value is not None and value>0]
        if positive and max(positive)/min(positive)>10:axes.set_yscale("log")
    fixed = {k: job["request"][k] for k in ("beta", "J0", "J1", "J2", "hx") if k != axis}
    axes.set_title(f'{job["request"]["Lx"]} x {job["request"]["Ly"]} periodic | '+", ".join(f"{k}={v:g}" for k, v in fixed.items()), fontsize=9)
    axes.legend(frameon=False, fontsize=7)
    figure.text(.02, .015, "Bars: 1 SE (tau: asymptotic estimate). Crosses require review; missing bars mean uncertainty unresolved.", fontsize=6)
    output = io.BytesIO()
    figure.savefig(output, format=format, dpi=300, metadata={"Creator": "SSE H2SQ benchmark"} if format != "png" else None)
    return output.getvalue()
