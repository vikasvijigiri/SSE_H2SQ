#!/usr/bin/env python3
"""
SSE parameter sweep generator
J0 fixed = 1
hx = T * tan(30°)
"""

from pathlib import Path
import shutil
import stat
import math

# ======================================================
# Global settings
# ======================================================

TOP_DIR = Path(__file__).resolve().parent

# Lattice sizes
L_LIST = [8, 12, 16, 20]

# Inverse temperatures
BETA_LIST = [4.0, 6.0, 8.0, 10.0]

# Couplings (J0 fixed = 1)
J1_LIST = [0.0, 0.1, 0.2, 0.3]
J2_LIST = [0.0, 0.05, 0.1]
J3_LIST = [0.0, 0.1]

# Monte Carlo parameters
N_THERM   = 100000
N_MEAS    = 300000


# Runtime files
JOBSCRIPT = TOP_DIR / "jobscript.sh"
MAKEFILE  = TOP_DIR / "Makefile"

TAN30 = math.tan(math.pi / 6.0)   # tan(30°)

# ======================================================
# Helpers
# ======================================================

def make_executable(path: Path):
    path.chmod(path.stat().st_mode | stat.S_IEXEC)


def write_input_file(
    path: Path,
    *,
    L,
    beta,
    J1, J2, J3,
    n_therm,
    n_measure,
):
    """Write SSE key-value input file"""

    T  = 1.0 / beta
    hx = T * TAN30

    with open(path, "w") as f:
        f.write(
            f"""# Lattice
            Lx        {L}
            Ly        {L}
            Lz        {L}

            # Temperature
            beta      {beta}

            # Couplings
            J0        1.0
            J1        {J1}
            J2        {J2}
            J3        {J3}
            hx        {hx:.8f}

            # Monte Carlo
            n_therm   {n_therm}
            n_measure {n_measure}
            """
        )


# ======================================================
# Main directory generation
# ======================================================

for beta in BETA_LIST:
    beta_dir = TOP_DIR / f".." / f"files" / f"beta_{beta:.2f}"
    beta_dir.mkdir(parents=True, exist_ok=True)

    for L in L_LIST:
        L_dir = beta_dir / f"L_{L}"
        L_dir.mkdir(exist_ok=True)

        for J1 in J1_LIST:
            for J2 in J2_LIST:
                for J3 in J3_LIST:

                    run_dir = L_dir / f"J1_{J1:.2f}_J2_{J2:.2f}_J3_{J3:.2f}"
                    run_dir.mkdir(exist_ok=True)

                    # Input file
                    write_input_file(
                        run_dir / "input_params.in",
                        L=L,
                        beta=beta,
                        J1=J1,
                        J2=J2,
                        J3=J3,
                        n_therm=N_THERM,
                        n_measure=N_MEAS
                    )

                    # Copy runtime files
                    shutil.copy(JOBSCRIPT, run_dir)
                    (run_dir / "Makefile").write_text(
                        MAKEFILE.read_text().replace("SDIR := ../src", "SDIR := ../../../../src")
                    )
                    make_executable(run_dir / JOBSCRIPT.name)

print("✔ SSE input files generated (J0 fixed, hx = T tan(30°))")
