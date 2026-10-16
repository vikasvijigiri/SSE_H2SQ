from pathlib import Path
import sys
BACKEND = Path(__file__).resolve().parents[2] / "backend"
sys.path.insert(0, str(BACKEND))
import sys
import os
import numpy as np

# Add ED analysis to path
sys.path.append(str(BACKEND / "Analysis" / "ED_analysis"))
from exact_diagonalzation_H2SQ import generate_black_plaquette_indices_4_4, generate_plaquette_neighbor_map, H2SQ_hamiltonian, xMag, zMag_square
import scipy.sparse.linalg as sla

Lx = Ly = 4
Nl = 1
beta = 1.0
J0 = 1.0
J1 = 0.3
J2 = 0.0
hx = 0.3

bplaqs = generate_black_plaquette_indices_4_4(Lx, Ly)
neighbor_map = generate_plaquette_neighbor_map(Lx, Ly)
Hsp = H2SQ_hamiltonian(Lx*Ly*Nl, J0, J1, J2, hx, bplaqs, neighbor_map)

W, V = sla.eigsh(Hsp, k=2, which='SA')

mx = xMag(Lx, Ly, Nl, beta, W, V)
m2, m4 = zMag_square(Lx, Ly, Nl, beta, W, V)

print(f"ED: Mx = {mx}")
print(f"ED: m2 = {m2}")
