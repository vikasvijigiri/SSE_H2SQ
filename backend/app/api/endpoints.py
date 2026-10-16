from fastapi import APIRouter, HTTPException
from fastapi.responses import PlainTextResponse
from app.services.qmc_results import import_legacy, jackknife_summary, result_text
from app.api.schemas import SimulationParams, SimulationResult
from app.services.simulation import run_dummy_simulation
from app.services import benchmark
from fastapi.responses import Response, FileResponse, JSONResponse
from starlette.concurrency import run_in_threadpool
import json
import asyncio

_figure_lock = asyncio.Lock()


router = APIRouter()


@router.get("/benchmarks")
async def benchmarks():
    return await run_in_threadpool(benchmark.list_jobs)


@router.post("/benchmarks", status_code=202)
async def start_benchmark(request: benchmark.BenchmarkRequest):
    try:
        return await benchmark.start_job(request)
    except RuntimeError as error:
        raise HTTPException(409, str(error)) from error


@router.post("/benchmarks/compare", status_code=201)
async def compare_benchmarks(identities: list[str]):
    try:
        return JSONResponse(await run_in_threadpool(benchmark.combine_jobs, identities), status_code=201)
    except FileNotFoundError as error:
        raise HTTPException(404, "Source experiment not found") from error
    except RuntimeError as error:
        raise HTTPException(409, str(error)) from error


@router.get("/benchmarks/{identity}")
async def benchmark_result(identity: str):
    # Files already contain JSON-native values. Bypass FastAPI's recursive
    # encoder for large multi-schedule ACF payloads polled by the dashboard.
    return JSONResponse(await _benchmark_data(identity))


async def _benchmark_data(identity: str):
    try:
        return await run_in_threadpool(benchmark.get_job, identity)
    except FileNotFoundError as error:
        raise HTTPException(404, "Benchmark not found") from error


@router.post("/benchmarks/{identity}/cancel")
async def cancel_benchmark(identity: str):
    return {"cancel_requested": benchmark.cancel_job(identity)}


@router.post("/benchmarks/{identity}/resume", status_code=202)
async def resume_benchmark(identity: str):
    try:
        return await benchmark.resume_job(identity)
    except FileNotFoundError as error:
        raise HTTPException(404, "Benchmark not found") from error
    except RuntimeError as error:
        raise HTTPException(409, str(error)) from error


@router.get("/benchmarks/{identity}/export")
async def export_benchmark(identity: str, format: str = "csv", axis: str = "beta", observable: str = "E", quantity: str = "mean"):
    job = await _benchmark_data(identity)
    if format == "csv":
        content, media = benchmark.csv_result(job), "text/csv"
    elif format == "json":
        content, media = json.dumps(job, indent=2, allow_nan=False), "application/json"
    elif format in {"pdf", "svg", "png"} and axis in job["request"]["axes"] and observable in benchmark.OBSERVABLES and quantity in {"mean", "difference", "tau_int", "ess_per_second"}:
        async with _figure_lock:
            content = await run_in_threadpool(benchmark.figure_result, job, axis, observable, format, quantity)
        media = {"pdf": "application/pdf", "svg": "image/svg+xml", "png": "image/png"}[format]
    else:
        raise HTTPException(422, "Choose csv, json, pdf, svg or png and a valid axis/observable")
    return Response(content, media_type=media, headers={"Content-Disposition": f'attachment; filename="benchmark-{identity[:8]}-{axis}-{observable}-{quantity}.{format}"'})


@router.get("/benchmarks/{identity}/traces/{point_id}/{method}/{chain}")
async def benchmark_trace(identity: str, point_id: str, method: str, chain: int):
    job = await _benchmark_data(identity)
    if method not in benchmark.METHODS or not 0 <= chain < job["request"]["chains"] or point_id not in {p["point_id"] for p in job["points"]}:
        raise HTTPException(404, "Trace not found")
    summary = next(p for p in job["points"] if p["point_id"] == point_id)["methods"].get(method, {})
    source_identity = summary.get("source_experiment", identity)
    source_point = summary.get("source_point", point_id)
    if source_identity != identity:
        source = await _benchmark_data(source_identity)
        if source_point not in {p["point_id"] for p in source["points"]}:
            raise HTTPException(404, "Source trace not found")
    path = benchmark.STORE / source_identity / source_point / method / str(chain) / "trace.npz"
    if not path.is_file():
        raise HTTPException(404, "This trace is not ready yet")
    return FileResponse(path, filename=f"qmc-{point_id}-{method}-chain{chain}.npz")


@router.post("/run-ed", response_model=SimulationResult)
async def run_ed(params: SimulationParams):
    try:
        result_data = await run_dummy_simulation(params, "ED")
    except ValueError as error:
        raise HTTPException(422, str(error)) from error
    return SimulationResult(
        status="success",
        message="ED Simulation completed",
        data=result_data
    )


@router.post("/run-qmc", response_model=SimulationResult)
async def run_qmc(params: SimulationParams):
    try:
        result_data = await run_dummy_simulation(params, "QMC")
    except ValueError as error:
        raise HTTPException(422, str(error)) from error
    return SimulationResult(
        status="success",
        message="QMC Simulation completed",
        data=result_data
    )

@router.post("/cancel/{client_id}")
async def cancel_simulation(client_id: str):
    from app.services.simulation import cancel_job
    found = cancel_job(client_id)
    return {"status": "cancelled" if found else "not_found", "client_id": client_id}


@router.get("/history")
async def get_history():
    from app.services.simulation import load_ed_cache, load_qmc_cache

    ed_cache = load_ed_cache()
    qmc_cache = load_qmc_cache()

    def _is_canonical_key(k: str) -> bool:
        return k.count("_") == 8  # 9 components → 8 underscores

    history_list = []
    all_keys = set(qmc_cache.keys()) | set(ed_cache.keys())

    for key in all_keys:
        if not _is_canonical_key(key):
            continue
            
        qmc_data = qmc_cache.get(key)
        ed_data = ed_cache.get(key)
        import_legacy(key, qmc_data)
        
        ts_qmc = qmc_data.get("timestamp") if qmc_data else None
        ts_ed = ed_data.get("timestamp") if ed_data else None
        
        ts = max(ts_qmc or "2000-01-01T00:00:00", ts_ed or "2000-01-01T00:00:00")
        
        history_list.append({
            "key": key,
            "ed": ed_data,
            "qmc": qmc_data,
            "qmc_jackknife": jackknife_summary(key),
            "timestamp": ts,
        })

    # Sort strictly by timestamp
    history_list.sort(key=lambda x: x["timestamp"], reverse=True)

    # Limit to last 5 runs
    history_list = history_list[:5]

    return {"status": "success", "data": history_list}


@router.get("/qmc-runs/{run_id}/result", response_class=PlainTextResponse)
async def get_qmc_result(run_id: str):
    text = result_text(run_id)
    if text is None:
        raise HTTPException(status_code=404, detail="Raw output is unavailable for this run.")
    return PlainTextResponse(text, headers={"Content-Disposition": f'inline; filename="SSE_data.txt"', "X-Content-Type-Options": "nosniff"})
