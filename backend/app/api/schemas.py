from typing import Optional
from pydantic import BaseModel, ConfigDict, Field, model_validator


class SimulationParams(BaseModel):
    model_config = ConfigDict(allow_inf_nan=False)
    # Common Parameters
    Lx: int = 4
    Ly: int = 4
    beta: float = Field(default=1.0, gt=0)
    J0: float = 1.0
    J1: float = 0.3
    J2: float = 0.025
    client_id: Optional[str] = None

    # ED Specific Parameters
    Nl: int = Field(default=1, ge=1)
    h: float = 0.0
    n: int = 250  # n_basis_vecs
    sparse: bool = True

    # QMC Specific Parameters
    Lz: int = Field(default=1, ge=1)
    J3: float = 0.0
    hx: float = 0.0
    n_therm: int = Field(default=500000, ge=0)
    n_measure: int = Field(default=500000, ge=1)
    debug_mode: bool = False

    @model_validator(mode="after")
    def validate_checkerboard(self):
        if self.Lx < 2 or self.Ly < 2 or self.Lx % 2 or self.Ly % 2:
            raise ValueError("Periodic checkerboard dimensions must be even and at least two.")
        return self



class SimulationResult(BaseModel):
    status: str
    message: str
    data: Optional[list] = None
