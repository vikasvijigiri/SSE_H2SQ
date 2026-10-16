from typing import Optional
from pydantic import BaseModel


class SimulationParams(BaseModel):
    # Common Parameters
    Lx: int = 4
    Ly: int = 4
    beta: float = 1.0
    J0: float = 1.0
    J1: float = 0.3
    J2: float = 0.025
    client_id: Optional[str] = None

    # ED Specific Parameters
    Nl: int = 1
    h: float = 0.0
    n: int = 250  # n_basis_vecs
    sparse: bool = True

    # QMC Specific Parameters
    Lz: int = 1
    J3: float = 0.0
    hx: float = 0.0
    n_therm: int = 500000
    n_measure: int = 500000
    debug_mode: bool = False



class SimulationResult(BaseModel):
    status: str
    message: str
    data: Optional[list] = None
