import asyncio
from typing import Optional
import subprocess
import json
import os
import sys
import shutil
import tempfile
import datetime
import pandas as pd

from app.api.schemas import SimulationParams

# Get exact path of exact_diagonalzation_H2SQ.py
# The directory of this file is backend/app/services/
# The Analysis folder is at backend/../Analysis
current_dir = os.path.dirname(os.path.abspath(__file__))
analysis_dir = os.path.abspath(os.path.join(current_dir, "..", "..", "Analysis", "ED_analysis"))
if analysis_dir not in sys.path:
    sys.path.append(analysis_dir)

try:
    from exact_diagonalzation_H2SQ import (
        estimate_all_observables, generate_black_plaquette_indices,
        H2SQ_hamiltonian, action, black_plaquette_neighbors,
        generate_plaquette_neighbor_map,
    )
except ImportError:
    from Analysis.ED_analysis.exact_diagonalzation_H2SQ import (
        estimate_all_observables, generate_black_plaquette_indices,
        H2SQ_hamiltonian, action, black_plaquette_neighbors,
        generate_plaquette_neighbor_map,
    )
# Registry of running QMC subprocesses keyed by client_id.
# Used by the /cancel endpoint to kill a job mid-run.
_active_processes: dict[str, subprocess.Popen] = {}

def cancel_job(client_id: str) -> bool:
    """Kill the QMC subprocess for the given client_id.  Returns True if found."""
    proc = _active_processes.get(client_id)
    if proc is None:
        return False
    try:
        proc.kill()   # kill() is cross-platform and more reliable than terminate()
    except Exception:
        pass
    return True

CACHE_DIR = os.path.abspath(os.path.join(current_dir, "..", "..", "files"))
CACHE_FILE = os.path.join(CACHE_DIR, "ed_cache.json")


def load_ed_cache() -> dict:
    if not os.path.exists(CACHE_FILE):
        return {}
    try:
        with open(CACHE_FILE, "r") as f:
            return json.load(f)
    except Exception:
        return {}


def save_ed_cache(cache: dict):
    os.makedirs(CACHE_DIR, exist_ok=True)
    try:
        with open(CACHE_FILE, "w") as f:
            json.dump(cache, f, indent=2)
    except Exception:
        pass


def get_ed_cache_key(params: SimulationParams) -> str:
    # Must match get_qmc_cache_key format (9 components) so history merging works.
    # ED never has J3; we hardcode 0.000000 in that slot.
    # Nl (ED) == Lz (QMC),  h (ED) == hx (QMC).
    return (
        f"{params.Lx}_{params.Ly}_{params.Nl}"
        f"_{params.J0:.6f}_{params.J1:.6f}_{params.J2:.6f}"
        f"_0.000000_{params.h:.6f}_{params.beta:.6f}"
    )


QMC_CACHE_FILE = os.path.join(CACHE_DIR, "qmc_cache.json")

def load_qmc_cache() -> dict:
    if not os.path.exists(QMC_CACHE_FILE):
        return {}
    try:
        with open(QMC_CACHE_FILE, "r") as f:
            return json.load(f)
    except Exception:
        return {}

def save_qmc_cache(cache: dict):
    os.makedirs(CACHE_DIR, exist_ok=True)
    try:
        with open(QMC_CACHE_FILE, "w") as f:
            json.dump(cache, f, indent=2)
    except Exception:
        pass

def get_qmc_cache_key(params: SimulationParams) -> str:
    # 9-component canonical key: Lx_Ly_Lz_J0_J1_J2_J3_hx_beta
    # Matches get_ed_cache_key when J3=0 and hx==h, Lz==Nl.
    return (
        f"{params.Lx}_{params.Ly}_{params.Lz}"
        f"_{params.J0:.6f}_{params.J1:.6f}_{params.J2:.6f}"
        f"_{params.J3:.6f}_{params.hx:.6f}_{params.beta:.6f}"
    )


def get_equivalent_ed_params(qmc_params: SimulationParams) -> SimulationParams:
    return SimulationParams(
        Lx=qmc_params.Lx,
        Ly=qmc_params.Ly,
        beta=qmc_params.beta,
        J0=qmc_params.J0,
        J1=qmc_params.J1,
        J2=qmc_params.J2,
        client_id=qmc_params.client_id,
        # ED specific
        Nl=qmc_params.Lz,
        h=qmc_params.hx,
        n=250,
        sparse=True
    )


def get_equivalent_qmc_params(ed_params: SimulationParams) -> SimulationParams:
    return SimulationParams(
        Lx=ed_params.Lx,
        Ly=ed_params.Ly,
        beta=ed_params.beta,
        J0=ed_params.J0,
        J1=ed_params.J1,
        J2=ed_params.J2,
        client_id=ed_params.client_id,
        # QMC specific
        Lz=ed_params.Nl,
        hx=ed_params.h,
        J3=ed_params.J3,
        n_therm=ed_params.n_therm,
        n_measure=ed_params.n_measure
    )


async def run_ed_simulation(params: SimulationParams):
    """Runs the actual Exact Diagonalization simulation in Python"""
    loop = asyncio.get_running_loop()
    spins = params.Lx * params.Ly * params.Nl
    client_id = params.client_id

    if spins > 16:
        raise ValueError(
            f"ED simulation is limited to a maximum of 16 spins to prevent memory overflow "
            f"(Current spins: {spins}). Please use QMC for larger system sizes."
        )

    def _send_time_est(phase: str, eta: float):
        """Send a time_estimate WebSocket message from the worker thread."""
        if not client_id:
            return
        
        async def _do_send():
            from app.main import active_connections
            # Wait up to 2 seconds for the websocket connection to be established
            for _ in range(20):
                if client_id in active_connections:
                    try:
                        await active_connections[client_id].send_text(json.dumps({
                            "type": "time_estimate",
                            "phase": phase,
                            "eta": round(eta, 2),
                        }))
                    except Exception:
                        pass
                    return
                await asyncio.sleep(0.1)
                
        asyncio.run_coroutine_threadsafe(_do_send(), loop)

    def _run():
        import time
        import itertools as _it
        import numpy as _np
        from scipy.linalg import eigh as _eigh
        from scipy.sparse.linalg import eigsh as _eigsh
        from exact_diagonalzation_H2SQ import basis_no

        bplaqs = generate_black_plaquette_indices(params.Lx, params.Ly)
        neighbor_map = generate_plaquette_neighbor_map(params.Lx, params.Ly)
        sparse = params.sparse
        n_vals = params.n
        dim = 2 ** spins

        if sparse:
            if dim <= 2:
                sparse = False
            else:
                n_vals = min(n_vals, dim - 2)
                if n_vals < 1:
                    n_vals = 1

        # ---- Phase 1: calibrate Hamiltonian build time ----
        # Replicate the FULL inner loop of H2SQ_hamiltonian (action + basis_no
        # + list appends) so the per-state cost matches exactly.
        try:
            CAL_N = min(500, dim)
            _data_cal, _rows_cal, _cols_cal = [], [], []
            t0 = time.perf_counter()
            for idx, ii in enumerate(_it.islice(_it.product([-1, 1], repeat=spins), CAL_N)):
                if idx % 100 == 0:
                    import time
                    time.sleep(0.001)
                for i, bplaq in enumerate(bplaqs):
                    nbh1 = bplaqs[neighbor_map[i][0]]
                    nbh2 = bplaqs[neighbor_map[i][1]]
                    state_list, coef_list = action(
                        ii, params.J0, params.J1, params.J2,
                        params.h, bplaq, nbh1, nbh2)
                    for jc, jj in enumerate(state_list):
                        _data_cal.append(coef_list[jc])
                        _rows_cal.append(basis_no(ii, spins))
                        _cols_cal.append(basis_no(jj, spins))
            cal_secs = time.perf_counter() - t0
            h_eta = (cal_secs / CAL_N) * max(0, dim - CAL_N)
        except Exception:
            h_eta = 0.0
        _send_time_est("hamiltonian", h_eta)

        # ---- Full Hamiltonian build ----
        Hsp = H2SQ_hamiltonian(spins, params.J0, params.J1, params.J2,
                               params.h, bplaqs, neighbor_map)

        # ---- Phase 2: estimate diagonalization time via trial matvecs ----
        # Lanczos (eigsh) costs ≈ n_vals * 3 matvecs; each matvec costs the
        # same as one Hsp.dot(v). Time 10 trial matvecs and scale.
        t_mv = 0.0
        try:
            v_trial = _np.random.randn(dim)
            N_MV = 10
            t_mv0 = time.perf_counter()
            for _ in range(N_MV):
                Hsp.dot(v_trial)
            t_mv = (time.perf_counter() - t_mv0) / N_MV
            if sparse:
                # Lanczos typically needs ~3–5× n_vals matvec iterations
                # Plus O(n_vals^2 * dim) FLOPs for orthogonalization
                ortho_time = (n_vals**2 * dim) / 1e8  # rough guesstimate of python/numpy speed
                diag_eta = (t_mv * n_vals * 4) + ortho_time
            else:
                # Dense eigh: O(dim^3); rough calibration from matvec cost
                diag_eta = t_mv * dim * 10
        except Exception:
            diag_eta = 0.0
        _send_time_est("diagonalization", diag_eta)

        if not sparse or dim <= 256:
            Wsp, Vsp = _eigh(Hsp.toarray())
        else:
            try:
                # Try default ncv first
                Wsp, Vsp = _eigsh(Hsp, k=n_vals, which="SA", return_eigenvectors=True)
            except Exception as e:
                # Catch ArpackError or other convergence issues
                if dim <= 2048:
                    # Fallback to dense for moderately small systems
                    Wsp, Vsp = _eigh(Hsp.toarray())
                    Wsp, Vsp = Wsp[:n_vals], Vsp[:, :n_vals]  # type: ignore
                else:
                    # Try again with a explicitly larger ncv and maxiter
                    try:
                        ncv_val = min(dim - 1, max(2 * n_vals + 1, 60))
                        Wsp, Vsp = _eigsh(Hsp, k=n_vals, which="SA", return_eigenvectors=True, ncv=ncv_val, maxiter=100000)
                    except Exception as e2:
                        raise RuntimeError(f"ED Diagonalization failed to converge: {str(e2)}")

        # ---- Phase 3: estimate observable calculation time ----
        try:
            # Observables are now fully vectorized and take ~O(dim) operations
            obs_eta = 0.5  # Fixed minimal ETA for fully vectorized observables
        except Exception:
            obs_eta = 0.0
        _send_time_est("observables", obs_eta)

        df = estimate_all_observables(
            params.Lx, params.Ly, params.Nl,
            params.J0, params.J1, params.J2,
            params.h, Wsp, Vsp, params.beta
        )
        return df.to_dict(orient="records")

    return await loop.run_in_executor(None, _run)


async def run_qmc_simulation(params: SimulationParams, client_id: Optional[str]):
    """Runs the compiled QMC C++ simulation executable"""
    build_dir = os.path.abspath(os.path.join(current_dir, "..", "..", "build"))
    debug_mode = getattr(params, 'debug_mode', False)
    if debug_mode:
        exe_name = "sse_sim_debug.exe" if os.name == "nt" else "sse_sim_debug"
        params = params.model_copy(update={"n_therm": 1, "n_measure": 0})
    else:
        exe_name = "sse_sim_test.exe" if os.name == "nt" else "sse_sim_test"
    exe_path = os.path.join(build_dir, exe_name)

    loop = asyncio.get_running_loop()

    async def _ws_send(msg: dict):
        """Send a WebSocket message from the async context, waiting briefly for the connection."""
        if not client_id:
            return
        from app.main import active_connections
        for _ in range(20):
            if client_id in active_connections:
                try:
                    await active_connections[client_id].send_text(json.dumps(msg))
                except Exception:
                    pass
                return
            await asyncio.sleep(0.1)

    if not os.path.exists(exe_path):
        await _ws_send({"type": "time_estimate", "phase": "compile", "eta": 60.0})
        try:
            def _run_make():
                make_target = "debugbin" if debug_mode else "all"
                return subprocess.run(
                    ["make", make_target], cwd=build_dir, capture_output=True, timeout=300
                )
            result = await loop.run_in_executor(None, _run_make)
            if result.returncode != 0:
                raise RuntimeError(
                    f"C++ Compilation failed:\n{result.stderr.decode(errors='replace')}"
                )
        except subprocess.TimeoutExpired:
            raise RuntimeError("C++ Compilation timed out after 300 seconds")
        await _ws_send({"type": "progress", "phase": "compile", "current": 1, "total": 1, "eta": 0.0})

    if not os.path.exists(exe_path):
        raise FileNotFoundError(f"QMC executable not found: {exe_path}")

    # Create a unique temp directory for this run
    workspace_root = os.path.abspath(os.path.join(current_dir, "..", ".."))
    temp_root = os.path.join(workspace_root, "files", "temp_runs")
    os.makedirs(temp_root, exist_ok=True)
    temp_dir = tempfile.mkdtemp(prefix=f"run_{client_id or 'api'}_", dir=temp_root)

    # ── Copy exe into temp dir so the build-dir binary is NEVER locked ──
    # On Windows, a running .exe is locked by the OS. By running a per-job
    # copy, we free the build-dir exe for future recompilation at any time.
    import shutil
    job_exe = os.path.join(temp_dir, exe_name)
    shutil.copy2(exe_path, job_exe)
    exe_path = job_exe  # use the per-job copy from here on

    try:
        # Write input file
        input_file_path = os.path.join(temp_dir, "input_params.in")
        with open(input_file_path, "w") as f:
            f.write(f"""# Lattice
Lx        {params.Lx}
Ly        {params.Ly}
Lz        {params.Lz}

# Temperature
beta      {params.beta}

# Couplings
J0        {params.J0}
J1        {params.J1}
J2        {params.J2}
J3        {params.J3}
hx        {params.hx}

# Monte Carlo
n_therm   {params.n_therm}
n_measure {params.n_measure}
""")
        
        # Run C++ QMC binary in a thread pool executor to support both Selector & Proactor loops on Windows
        def _execute_qmc():
            from app.main import active_connections

            proc = subprocess.Popen(
                [exe_path, "input_params.in"],
                cwd=temp_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                encoding='utf-8',
                errors='replace',
                bufsize=1
            )

            # Register so the /cancel endpoint can kill us.
            if client_id:
                _active_processes[client_id] = proc

            # Accumulate SSE algorithm check results from [SSE-CHECK] stdout lines.
            check_counts: dict = {}   # {name: {"total": int, "passed": int}}
            last_energy_inst: list = [None]  # list wrapper to allow mutation in loop
            summary_total: list = [0]
            summary_failed: list = [0]

            # Read stdout line by line; C++ emits [PROGRESS] and [SSE-CHECK] lines.
            try:
                if proc.stdout:
                    while True:
                        line = proc.stdout.readline()
                        if not line:
                            break
                        # Stream every raw stdout line to the terminal viewer in the UI
                        stripped_line = line.rstrip('\n').rstrip('\r')
                        if stripped_line and client_id and client_id in active_connections:
                            asyncio.run_coroutine_threadsafe(
                                active_connections[client_id].send_text(json.dumps({
                                    "type": "stdout",
                                    "line": stripped_line,
                                })),
                                loop
                            )
                        if line.startswith("[PROGRESS]"):
                            try:
                                # [PROGRESS] therm 250000 500000 12.50
                                parts = line.strip().split()
                                if len(parts) >= 5:
                                    phase   = parts[1]
                                    current = int(parts[2])
                                    total   = int(parts[3])
                                    eta_val = float(parts[4])
                                    if client_id and client_id in active_connections:
                                        ws = active_connections[client_id]
                                        asyncio.run_coroutine_threadsafe(
                                            ws.send_text(json.dumps({
                                                "type": "progress",
                                                "phase": phase,
                                                "current": current,
                                                "total": total,
                                                "eta": eta_val
                                            })),
                                            loop
                                        )
                            except Exception as e:
                                print("Error parsing progress:", e)
                        elif line.startswith("[SSE-CHECK]"):
                            try:
                                # [SSE-CHECK] PASS n_ops step=50000 phase=therm
                                # [SSE-CHECK] FAIL n_ops step=50000 phase=therm
                                # [SSE-CHECK] INFO energy step=5000 phase=measure n_ops=47 E_inst=-0.5234
                                parts = line.strip().split()
                                if len(parts) >= 3:
                                    status = parts[1]   # PASS | FAIL | INFO
                                    name   = parts[2]   # n_ops | weights | tau_period | energy

                                    if status in ("PASS", "FAIL"):
                                        if name not in check_counts:
                                            check_counts[name] = {"total": 0, "passed": 0}
                                        check_counts[name]["total"] += 1
                                        passed = (status == "PASS")
                                        if passed:
                                            check_counts[name]["passed"] += 1

                                        # Extract step/phase for the WS message
                                        step_val, phase_val = None, None
                                        for kv in parts[3:]:
                                            if kv.startswith("step="):
                                                try: step_val = int(kv.split("=", 1)[1])
                                                except: pass
                                            elif kv.startswith("phase="):
                                                phase_val = kv.split("=", 1)[1]

                                        if client_id and client_id in active_connections:
                                            ws = active_connections[client_id]
                                            asyncio.run_coroutine_threadsafe(
                                                ws.send_text(json.dumps({
                                                    "type": "check_result",
                                                    "name": name,
                                                    "passed": passed,
                                                    "step": step_val,
                                                    "phase": phase_val,
                                                })),
                                                loop
                                            )

                                    elif status == "INFO" and name == "energy":
                                        for kv in parts[3:]:
                                            if kv.startswith("E_inst="):
                                                try: last_energy_inst[0] = float(kv.split("=", 1)[1])
                                                except: pass

                            except Exception as e:
                                print("Error parsing SSE-CHECK:", e)
                        elif line.startswith("[SSE-CHECKS-SUMMARY]"):
                            try:
                                # [SSE-CHECKS-SUMMARY] total=30 failed=0
                                for kv in line.strip().split()[1:]:
                                    if kv.startswith("total="):
                                        try: summary_total[0] = int(kv.split("=", 1)[1])
                                        except: pass
                                    elif kv.startswith("failed="):
                                        try: summary_failed[0] = int(kv.split("=", 1)[1])
                                        except: pass
                            except Exception as e:
                                print("Error parsing SSE-CHECKS-SUMMARY:", e)
            finally:
                proc.wait()
                if client_id:
                    _active_processes.pop(client_id, None)

            stderr_data = proc.stderr.read() if proc.stderr else ""

            # Build per-check summary dict
            checks_summary: dict = {}
            for cname, counts in check_counts.items():
                checks_summary[cname] = "PASS" if counts["passed"] == counts["total"] else "FAIL"
            total_c  = summary_total[0]  or sum(c["total"]           for c in check_counts.values())
            failed_c = summary_failed[0] or sum(c["total"] - c["passed"] for c in check_counts.values())
            if total_c > 0:
                checks_summary["_total"]      = total_c
                checks_summary["_failed"]     = failed_c
                checks_summary["_all_passed"] = (failed_c == 0)
            if last_energy_inst[0] is not None:
                checks_summary["_energy_inst"] = round(last_energy_inst[0], 6)

            return proc.returncode, stderr_data, checks_summary

        returncode, stderr_data, checks_summary = await loop.run_in_executor(None, _execute_qmc)
        
        if returncode != 0:
            raise RuntimeError(f"QMC C++ simulation failed: {stderr_data}")
            
        # Parse result data from SSE_data.txt
        output_file_path = os.path.join(temp_dir, "SSE_data.txt")
        if not os.path.exists(output_file_path):
            raise FileNotFoundError("Simulation finished but output file SSE_data.txt was not created.")
            
        df = pd.read_csv(output_file_path, sep=r'\s+')
        records = df.to_dict(orient="records")
        if not records:
            raise ValueError("No observable records found in the output data file.")

        last_row = records[-1]

        # Parse spin configuration snapshots written by the C++ binary
        spin_configs = []
        spin_config_path = os.path.join(temp_dir, "spin_configs.json")
        if os.path.exists(spin_config_path):
            try:
                with open(spin_config_path, "r") as f:
                    spin_configs = json.load(f)
            except Exception:
                spin_configs = []

        mapped_data = [{
            "Lx": params.Lx,
            "Ly": params.Ly,
            "Lz": params.Lz,
            "beta": last_row.get("beta", params.beta),
            "J0": last_row.get("J0", params.J0),
            "J1": last_row.get("J1", params.J1),
            "J2": last_row.get("J2", params.J2),
            "J3": last_row.get("J3", params.J3),
            "hx": last_row.get("hx", params.hx),
            "n_therm": params.n_therm,
            "n_measure": params.n_measure,
            # Energy moments
            "enrg":   last_row.get("E",  0.0),
            "enrg2":  last_row.get("E2", 0.0),
            "enrg4":  last_row.get("E4", 0.0),
            "E_trans": last_row.get("E_trans", 0.0),
            # z-magnetization moments
            "mz":          last_row.get("mz",  0.0),
            "SMag_square": last_row.get("mz2", 0.0),
            "SMag_four":   last_row.get("mz4", 0.0),
            # x-magnetization moments
            "Mag_x": last_row.get("mx",  0.0),
            "mx2":   last_row.get("mx2", 0.0),
            "mx4":   last_row.get("mx4", 0.0),
            # Order parameter moments
            "OP":  last_row.get("OP",  0.0),
            "OP2": last_row.get("OP2", 0.0),
            "OP4": last_row.get("OP4", 0.0),
            # Ice-rule observable
            "ice":  last_row.get("ice",  0.0),
            "ice2": last_row.get("ice2", 0.0),
            # Superfluid stiffness (helicity modulus)
            "rho_x": last_row.get("rho_x", 0.0),
            "rho_y": last_row.get("rho_y", 0.0),
            "spin_configs": spin_configs,
            "checks_summary": checks_summary if checks_summary else None,
            "timestamp": datetime.datetime.now().isoformat()
        }]
        return mapped_data
        
    finally:
        try:
            shutil.rmtree(temp_dir)
        except Exception:
            pass

async def run_dummy_simulation(params: SimulationParams, sim_type: str):
    """Wrapper function to execute the real simulation and stream updates over WebSockets"""
    from app.main import active_connections
    client_id = params.client_id
    
    # 1. Prepare param objects for both
    if sim_type == "ED":
        ed_params = params
        qmc_params = get_equivalent_qmc_params(params)
    else:
        qmc_params = params
        ed_params = get_equivalent_ed_params(params)

    debug_mode = sim_type == "QMC" and getattr(qmc_params, 'debug_mode', False)

    # Ensure ONLY the primary simulation reports progress to the UI
    ed_client_id = client_id if sim_type == "ED" else None
    qmc_client_id = client_id if sim_type == "QMC" else None
    
    ed_params.client_id = ed_client_id
    qmc_params.client_id = qmc_client_id

    # 2. Check caches
    ed_cache_key = get_ed_cache_key(ed_params)
    ed_cache = load_ed_cache()
    ed_record = ed_cache.get(ed_cache_key)

    qmc_cache_key = get_qmc_cache_key(qmc_params)
    qmc_cache = load_qmc_cache()
    qmc_record = qmc_cache.get(qmc_cache_key)

    # STRICT MATCHING: ED does not support J3.
    # If the physical parameters include J3 != 0, ED cannot represent this state.
    if abs(qmc_params.J3) > 1e-6:
        ed_record = None

    # If the user explicitly launches QMC, do not use the cache because QMC is stochastic.
    if sim_type == "QMC":
        qmc_record = None

    # 3. Create tasks for missing records
    run_ed_task = None
    run_qmc_task = None
    
    spins_ed = ed_params.Lx * ed_params.Ly * ed_params.Nl
    
    if not ed_record and abs(qmc_params.J3) <= 1e-6 and not debug_mode:
        if sim_type == "ED" and spins_ed <= 16:
            run_ed_task = asyncio.create_task(run_ed_simulation(ed_params))
        elif sim_type == "QMC" and spins_ed <= 16:
            run_ed_task = asyncio.create_task(run_ed_simulation(ed_params))

    if not qmc_record:
        # QMC can run any size, handle background caching seamlessly
        run_qmc_task = asyncio.create_task(run_qmc_simulation(qmc_params, qmc_client_id))

    # 4. Wait for the PRIMARY requested task
    primary_task = run_ed_task if sim_type == "ED" else run_qmc_task
    
    if primary_task:
        try:
            res = await primary_task
            if sim_type == "ED" and res:
                ed_record = res[-1]
                import datetime
                ed_record["timestamp"] = datetime.datetime.now().isoformat()
                ed_record.setdefault("Lz", ed_record.get("Nl", ed_params.Nl))
                ed_record.setdefault("hx", ed_record.get("h", ed_params.h))
                ed_record.setdefault("J3", 0.0)
                c = load_ed_cache()
                c[ed_cache_key] = ed_record
                save_ed_cache(c)
            elif sim_type == "QMC" and res:
                qmc_record = res[-1]
                if not debug_mode:
                    import datetime
                    qmc_record["timestamp"] = datetime.datetime.now().isoformat()
                    c = load_qmc_cache()
                    c[qmc_cache_key] = qmc_record
                    save_qmc_cache(c)
        except Exception as e:
            raise e

    # 5. For secondary (background) task: wait only briefly (5s) then let it run in background.
    # This way QMC always returns promptly. ED result will be in cache for next QMC run.
    secondary_task = run_ed_task if sim_type == "QMC" else None
    if secondary_task:
        try:
            res2 = await asyncio.wait_for(secondary_task, timeout=5.0)
            if res2:
                import datetime
                ed_record = res2[-1]
                ed_record["timestamp"] = datetime.datetime.now().isoformat()
                ed_record.setdefault("Lz", ed_record.get("Nl", ed_params.Nl))
                ed_record.setdefault("hx", ed_record.get("h", ed_params.h))
                ed_record.setdefault("J3", 0.0)
                c = load_ed_cache()
                c[ed_cache_key] = ed_record
                save_ed_cache(c)
        except (asyncio.TimeoutError, Exception):
            pass  # ED is slow or failed — QMC result still returns immediately

    # Reload caches to get any new background results
    ed_cache = load_ed_cache()
    qmc_cache = load_qmc_cache()
    ed_record = ed_cache.get(ed_cache_key)
    if not debug_mode:
        qmc_record = qmc_cache.get(qmc_cache_key)

    results = [{"qmc": qmc_record, "ed": ed_record}]
        
    # Stream final 100% progress and complete event
    if client_id and client_id in active_connections:
        ws = active_connections[client_id]
        try:
            await ws.send_text(json.dumps({
                "type": "progress",
                "sim_type": sim_type,
                "progress": 100,
                "eta": 0.0
            }))
            await ws.send_text(json.dumps({
                "type": "complete",
                "sim_type": sim_type
            }))
        except Exception:
            pass
            
    return results
