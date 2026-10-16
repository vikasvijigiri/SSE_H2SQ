from fastapi import APIRouter, HTTPException
from fastapi.responses import PlainTextResponse
from app.services.qmc_results import import_legacy, jackknife_summary, result_text
from app.api.schemas import SimulationParams, SimulationResult
from app.services.simulation import run_dummy_simulation


router = APIRouter()


@router.post("/run-ed", response_model=SimulationResult)
async def run_ed(params: SimulationParams):
    result_data = await run_dummy_simulation(params, "ED")
    return SimulationResult(
        status="success",
        message="ED Simulation completed",
        data=result_data
    )


@router.post("/run-qmc", response_model=SimulationResult)
async def run_qmc(params: SimulationParams):
    result_data = await run_dummy_simulation(params, "QMC")
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
