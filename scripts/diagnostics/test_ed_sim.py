from pathlib import Path
import sys
BACKEND = Path(__file__).resolve().parents[2] / "backend"
sys.path.insert(0, str(BACKEND))
import traceback
import asyncio
from app.services.simulation import run_ed_simulation
from app.api.schemas import SimulationParams

params = SimulationParams(Lx=4, Ly=4, Nl=1, sparse=True)

async def main():
    try:
        df = await run_ed_simulation(params)
        print('SUCCESS')
    except Exception as e:
        traceback.print_exc()

asyncio.run(main())
