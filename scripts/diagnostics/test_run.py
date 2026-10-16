from pathlib import Path
import sys
BACKEND = Path(__file__).resolve().parents[2] / "backend"
sys.path.insert(0, str(BACKEND))
import asyncio
from app.api.schemas import SimulationParams
from app.services.simulation import run_dummy_simulation

async def test():
    params = SimulationParams(
        Lx=4, Ly=4, Nl=1, J0=1.0, J1=0.3, J2=0.025, h=0.0, beta=10.0, n=250, sparse=True
    )
    # The UI shows "Beta=10" in the screenshot.
    # Notice that `Beta=10` is used! For ED, beta=10 might be very cold.
    print("Testing ED...")
    res = await run_dummy_simulation(params, "ED")
    print(res)

if __name__ == "__main__":
    asyncio.run(test())
