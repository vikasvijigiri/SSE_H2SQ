from pathlib import Path
import sys
BACKEND = Path(__file__).resolve().parents[2] / "backend"
sys.path.insert(0, str(BACKEND))
import asyncio
from app.api.schemas import SimulationParams
from app.services.simulation import run_qmc_simulation

async def test():
    params = SimulationParams(
        Lx=4, Ly=4, Lz=1, J0=1.0, J1=0.3, J2=0.025, hx=0.0, beta=10.0, n_therm=1000, n_measure=1000
    )
    print("Testing QMC...")
    try:
        res = await run_qmc_simulation(params, None)
        print(res)
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    asyncio.run(test())
