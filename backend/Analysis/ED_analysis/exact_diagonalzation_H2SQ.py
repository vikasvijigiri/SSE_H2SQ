# ***********************************************************************************
#                              Exact diagonalization 
# Description: Spin-1/2 Bilayer model on square lattice
# Author: Vikas Vijigiri
# Date: October, 2023
# ***********************************************************************************
import numpy as np
from   scipy.sparse import csr_matrix
import scipy as sp
from scipy.special import logsumexp
from   scipy.sparse.linalg import eigsh
import itertools as it
import pandas as pd
import concurrent.futures
import os
from typing import List, Optional

H_BONDS = 4
Q_BONDS = 4
BQ_BONDS = 2



# ***********************************************************************************
#                               Hamiltonian part 
# ***********************************************************************************


class Plaquette:
    def __init__(self, indices):
        self.indices = indices

def mod(a, b):
    return (a % b + b) % b

def generate_lattice_indices(Lx, Ly):
    lattice_indices = [[0] * Ly for _ in range(Lx)]

    for i in range(Lx):
        for j in range(Ly):
            lattice_indices[i][j] = i * Ly + j

    return lattice_indices

def generate_black_plaquette_indices(Lx, Ly):
    lattice_indices = generate_lattice_indices(Lx, Ly)
    black_plaquette_indices = []

    for i in range(Lx):
        for j in range(Ly):
            if (i + j) % 2 == 0:
                plaquette = Plaquette([
                    lattice_indices[i][j],
                    lattice_indices[i][mod((j + 1), Ly)],
                    lattice_indices[mod((i + 1), Lx)][mod((j + 1), Ly)],
                    lattice_indices[mod((i + 1), Lx)][j]
                ])
                black_plaquette_indices.append(plaquette)

    return black_plaquette_indices

def generate_black_plaquette_indices_4_4(Lx, Ly):
    lattice_indices = generate_lattice_indices(Lx, Ly)
    black_plaquette_indices = []

#    for i in range(Lx - 1):
#        for j in range(Ly - 1):
#            if (i + j) % 2 == 0:
#                plaquette = Plaquette([
#                    lattice_indices[i][j],
#                    lattice_indices[i][j + 1],
#                    lattice_indices[i + 1][j + 1],
#                    lattice_indices[i + 1][j]
#                ])
#                black_plaquette_indices.append(plaquette)

    plaquette = Plaquette([
        0,
        1,
        4,
        5
    ])
    black_plaquette_indices.append(plaquette)

    plaquette = Plaquette([
        2,
        3,
        6,
        7
    ])
    black_plaquette_indices.append(plaquette)

    plaquette = Plaquette([
        5,
        6,
        9,
        10
    ])
    black_plaquette_indices.append(plaquette)

    plaquette = Plaquette([
        7,
        4,
        11,
        8
    ])
    black_plaquette_indices.append(plaquette)    

    plaquette = Plaquette([
        8,
        9,
        12,
        13
    ])
    black_plaquette_indices.append(plaquette)       

    plaquette = Plaquette([
        10,
        11,
        14,
        15
    ])
    black_plaquette_indices.append(plaquette)       


    plaquette = Plaquette([
        13,
        14,
        1,
        2
    ])
    black_plaquette_indices.append(plaquette)     

    plaquette = Plaquette([
        15,
        12,
        3,
        0
    ])
    black_plaquette_indices.append(plaquette)             

    return black_plaquette_indices

def generate_black_plaquette_indices_4_2(Lx, Ly):
    lattice_indices = generate_lattice_indices(Lx, Ly)
    black_plaquette_indices = []

#    for i in range(Lx - 1):
#        for j in range(Ly - 1):
#            if (i + j) % 2 == 0:
#                plaquette = Plaquette([
#                    lattice_indices[i][j],
#                    lattice_indices[i][j + 1],
#                    lattice_indices[i + 1][j + 1],
#                    lattice_indices[i + 1][j]
#                ])
#                black_plaquette_indices.append(plaquette)

    plaquette = Plaquette([
        0,
        1,
        4,
        5
    ])
    black_plaquette_indices.append(plaquette)

    plaquette = Plaquette([
        2,
        3,
        6,
        7
    ])
    black_plaquette_indices.append(plaquette)

    return black_plaquette_indices

def generate_black_plaquette_indices_4_2_2(Lx, Ly):
    lattice_indices = generate_lattice_indices(Lx, Ly)
    black_plaquette_indices = []

#    for i in range(Lx - 1):
#        for j in range(Ly - 1):
#            if (i + j) % 2 == 0:
#                plaquette = Plaquette([
#                    lattice_indices[i][j],
#                    lattice_indices[i][j + 1],
#                    lattice_indices[i + 1][j + 1],
#                    lattice_indices[i + 1][j]
#                ])
#                black_plaquette_indices.append(plaquette)

    plaquette = Plaquette([
        0,
        1,
        4,
        5
    ])
    black_plaquette_indices.append(plaquette)

    plaquette = Plaquette([
        2,
        3,
        6,
        7
    ])

    plaquette = Plaquette([
        8,
        9,
        12,
        13
    ])
    black_plaquette_indices.append(plaquette)

    return black_plaquette_indices


def generate_plaquette_neighbor_map(Lx, Ly):
    """
    Build a dict mapping each black plaquette index to its 4 diagonal
    neighbors [BL, BR, TL, TR] for any Lx x Ly lattice with PBC.
    """
    pos_to_idx = {}
    positions = []
    for x in range(Lx):
        for y in range(Ly):
            if (x + y) % 2 == 0:
                pos_to_idx[(x, y)] = len(positions)
                positions.append((x, y))
    neighbor_map = {}
    for p_idx, (x, y) in enumerate(positions):
        nbrs = []
        for dx, dy in [(-1, -1), (1, -1), (-1, 1), (1, 1)]:
            nx, ny = (x + dx) % Lx, (y + dy) % Ly
            nbrs.append(pos_to_idx[(nx, ny)])
        neighbor_map[p_idx] = nbrs
    return neighbor_map


def black_plaquette_neighbors(index):
    """
    Return the 4 neighbor plaquette indices for a given black plaquette
    index on a 4x4 lattice (hardcoded for backward compatibility).
    Use generate_plaquette_neighbor_map(Lx, Ly) for other lattice sizes.
    """
    neighbor_map = {
        0: [7, 6, 3, 2],  # bottom left, bottom right, top left, top right
        1: [6, 7, 2, 3],
        2: [0, 1, 4, 5],
        3: [1, 0, 5, 4],
        4: [3, 2, 7, 6],
        5: [2, 3, 6, 7],
        6: [4, 5, 0, 1],
        7: [5, 4, 1, 0]
    }
    return neighbor_map[index]



def print_indices(plaquettes):
    for plaquette in plaquettes:
        for index in plaquette.indices:
            print(index, end=" ")
        print()

## Example usage:
#Lx = 4
#Ly = 4

#black_plaquettes = generate_black_plaquette_indices(Lx, Ly)
#print("Black Plaquette Indices:")
#print_indices(black_plaquettes)



def basis_no(basis_state, L):
    return sum((2**j) * ((basis_state[j] + 1) // 2) for j in range(L))


def action(ii, J0, J1, J2, h, bplaq, bplaq_nbh1, bplaq_nbh2):
    i, j, k, l = bplaq.indices
    i1, j1, k1, l1 = bplaq_nbh1.indices
    i2, j2, k2, l2 = bplaq_nbh2.indices    
    
    state_list = []
    coef_list = []
    
    # Heisenberg AFM only model
    # diagonal
    de  = - J0 * (ii[i] * ii[k] * ii[j] * ii[l])
    # Physical diagonal pair interaction on cyclic plaquette sites.
    de += J1 * (ii[i]*ii[k] + ii[j]*ii[l])
    de +=   J2 * (ii[i1] + ii[l1] - ii[j1] - ii[k1]) / 4.0 * (ii[i] + ii[l] - ii[j] - ii[k]) / 4.0
    de +=   J2 * (ii[i1] + ii[j1] - ii[k1] - ii[l1]) / 4.0 * (ii[i] + ii[j] - ii[k] - ii[l]) / 4.0
    de +=   J2 * (ii[i2] + ii[l2] - ii[j2] - ii[k2]) / 4.0 * (ii[i] + ii[l] - ii[j] - ii[k]) / 4.0
    de +=   J2 * (ii[i2] + ii[j2] - ii[k2] - ii[l2]) / 4.0 * (ii[i] + ii[j] - ii[k] - ii[l]) / 4.0    
    if abs(de) > 1e-10:
      coef_list.append(de)
      state_list.append(ii)

    # off-diagonal
    for s in [i, j, k, l]:
        jj = list(ii)
        jj[s] = -ii[s] # flip site i
        jj = tuple(jj)
        coef_list.append(-0.5 * h)     # -0.5 * h because each site is in 2 black plaquettes
        state_list.append(jj)

    return state_list, coef_list   

# Spin-1/2
def H2SQ_hamiltonian(L, J0, J1, J2, h, bplaqs, neighbor_map=None):

    data = []
    rows = []
    cols = []

    if neighbor_map is None:
        neighbor_map = black_plaquette_neighbors

    for idx, ii in enumerate(it.product([-1, 1], repeat=L)):
        if idx % 1024 == 0:
            import time
            time.sleep(0.001)
        # Ising (J0 + J1 + J2) + zeeman (h) interaction interactions
        for i, bplaq in enumerate(bplaqs):
            nbrs = neighbor_map(i) if callable(neighbor_map) else neighbor_map[i]
            bplaq_nbh1 = bplaqs[nbrs[0]]
            bplaq_nbh2 = bplaqs[nbrs[1]]
            state_list, coef_list = action(ii, J0, J1, J2, h, bplaq, bplaq_nbh1, bplaq_nbh2)
            for jc, jj in enumerate(state_list):
                data.append(coef_list[jc])
                rows.append(basis_no(ii, L))
                cols.append(basis_no(jj, L)) 

    Hsp = csr_matrix((data, (rows, cols)), shape=(2**L, 2**L), dtype=np.float64)     
    return Hsp 



def H2SQ_ED(Lx, Ly, Nl, J0, J1, J2, h, Bplqs, n, sparse=False):
    #if n < 30: sparse = True
    Hsp = H2SQ_hamiltonian(Lx*Ly*Nl, J0, J1, J2, h, Bplqs)
    print("Hamiltonian Built -> success!")
    if not sparse: Wsp, Vsp = sp.linalg.eigh(Hsp.toarray())
    else: Wsp, Vsp = eigsh(Hsp, k=n, which='SA',return_eigenvectors=True)
    print("Diagonalization -> success!")
    return Wsp, Vsp


# ***********************************************************************************
#                               Measurements 
# ***********************************************************************************

class Site:
    def __init__(self, x, y):
        self.x = x
        self.y = y

def calculate_position_vectors(Lx, Ly):
    positions = []

    # Reference position for site 0
    ref_x, ref_y = 0, 0

    # Calculate position vectors
    for i in range(Lx * Ly):
        rel_x = i % Lx
        rel_y = i // Lx
        x = rel_x - ref_x
        y = rel_y - ref_y

        # Apply periodic boundary conditions along x
        if x > Lx / 2:
            x -= Lx
        if x < -Lx / 2:
            x += Lx

        # Apply periodic boundary conditions along y
        if y > Ly / 2:
            y -= Ly
        if y < -Ly / 2:
            y += Ly

        positions.append(Site(x, y))

    return positions


# def zMag_square(Lx, Ly, Nl, beta, basis_list, W, V):
#     smag_square_L1 = smag_four_L1 = 0.
#
#     z = np.sum(np.exp(-np.longdouble(beta * W)))
#
#     for ii, w in enumerate(W):
#         state_vector = V[:,ii]
#         weight = np.exp(-np.longdouble(beta * w))
#         for ii, basis in enumerate(basis_list):
#             si = state_vector[ii] ** 2
#             mg1 = 0. # layer 1
#             for x in range(Lx*Ly):
#                 mg1 += 0.5 * basis[x] * (-1) ** ((x // Lx) + (x % Lx))
#
#             smag_square_L1 += ((mg1/Lx/Ly/Nl)**2) * si * weight/z
#             smag_four_L1 += ((mg1/Lx/Ly/Nl)**4) * si * weight/z
#
#     return smag_square_L1, smag_four_L1

def zMag_square(Lx, Ly, Nl, beta, W, V):
    Ns = Lx * Ly * Nl
    dim = 2**Ns

    # --- stable Boltzmann weights ---
    betaW = beta * W
    bmin = np.min(betaW)
    weights = np.exp(-(betaW - bmin))
    z = np.sum(weights)

    smag_square_L1 = 0.0
    smag_four_L1   = 0.0

    # --- precompute staggered magnetization for each basis ---
    stagger = np.zeros(dim)
    for s in range(dim):
        m = 0.0
        for x in range(Ns):
            spin = 1 if (s & (1 << x)) else -1
            m += 0.5 * spin * (-1) ** ((x // Lx) + (x % Lx))
        stagger[s] = m / Ns

    # --- main sum ---
    for eig_i, w in enumerate(W):
        state_vector = V[:, eig_i]
        weight = weights[eig_i] / z

        probs = state_vector**2   # |ψ|²
        smag_square_L1 += np.sum((stagger**2) * probs) * weight
        smag_four_L1   += np.sum((stagger**4) * probs) * weight

    return smag_square_L1, smag_four_L1



def xMag(Lx, Ly, Nl, beta, W, V):
    """
    Compute <Mx> for spin-1/2 system using ED
    """
    Ns = Lx * Ly * Nl
    dim = 2**Ns

    # ---- stable Boltzmann weights ----
    betaW = beta * W
    bmin = np.min(betaW)
    weights = np.exp(-(betaW - bmin))
    Z = np.sum(weights)

    Mx = 0.0
    indices = np.arange(dim)

    # ---- eigenstate loop ----
    for n, En in enumerate(W):
        psi = V[:, n]
        weight = weights[n] / Z

        mx_n = 0.0

        # ---- site loop ----
        for site in range(Ns):
            flipped_indices = indices ^ (1 << site)
            contrib = np.sum(psi * psi[flipped_indices])
            mx_n += contrib

        mx_n /= Ns
        Mx += weight * mx_n

    return Mx


# def zMag_square(Lx, Ly, Nl, beta, basis_list, W, V):
#     Ns = Lx * Ly

#     # energies per site
#     E = W / (Lx * Ly * Nl)

#     logZ = logsumexp(-beta * W)
#     weights = np.exp(-beta * W - logZ)

#     # precompute staggered magnetization
#     stagger = np.zeros(len(basis_list))
#     for bi, basis in enumerate(basis_list):
#         m = 0.0
#         for x in range(Ns):
#             m += 0.5 * basis[x] * (-1) ** ((x // Lx) + (x % Lx))
#         stagger[bi] = m / (Lx * Ly * Nl)

#     smag_square_L1 = 0.0
#     smag_four_L1   = 0.0

#     for eig_i in range(len(W)):
#         probs = V[:, eig_i]**2
#         smag_square_L1 += np.sum(stagger**2 * probs) * weights[eig_i]
#         smag_four_L1   += np.sum(stagger**4 * probs) * weights[eig_i]

#     return smag_square_L1, smag_four_L1


def estimate_all_observables(Lx, Ly, Nl, J0, J1, J2, h, Wsp, Vsp, beta):
    # # Rescale energies per site
    # E = Wsp / (Lx * Ly * Nl)

    # # Log-partition function (stable)
    # logZ = logsumexp(-beta * Wsp)

    # # Normalized Boltzmann weights (stable)
    # weights = np.exp(-beta * Wsp - logZ)

    # # Observables
    # enrg  = np.sum(E    * weights)
    # enrg2 = np.sum(E**2 * weights)
    # enrg4 = np.sum(E**4 * weights)

    # Rescale energies per site
    E = Wsp / (Lx * Ly * Nl)

    # Stable Boltzmann weights
    betaW = beta * Wsp
    bmin = np.min(betaW)

    weights = np.exp(-(betaW - bmin))   # all exponents ≤ 0
    z = np.sum(weights)

    # Observables
    enrg  = np.sum(E      * weights) / z
    enrg2 = np.sum(E**2   * weights) / z
    enrg4 = np.sum(E**4   * weights) / z




    basis_list = list(it.product([-1, 1], repeat=Lx*Ly*Nl))
    
    
    # Staggered Magnetization, Staggered Magnetization^2, Staggered Magnetization^4
    SMag_square, SMag_four = zMag_square(Lx, Ly, Nl, beta, Wsp, Vsp)
    mag_x = xMag(Lx, Ly, Nl, beta, Wsp, Vsp)


    # Printing the observables in one go using pandas Dataframe (With clear headers)
    headers = [ "Lx", "Ly", "Nl", "beta", "J0", "J1", "J2", "h", "n_basis_vecs",
                "enrg", "enrg2", "enrg4", 
                "SMag_square", "SMag_four", "Mag_x"
              ]

    values  = [ Lx, Ly, Nl, beta, J0, J1, J2, h, len(Wsp),
                enrg, enrg2, enrg4,
                SMag_square, SMag_four, mag_x 
              ]

    values_reshaped = [[val] for val in values]

    return pd.DataFrame([values], columns=headers)



def dump_df(df: pd.DataFrame, filename: str = "ED_data.csv", append: bool = True):
    """
    Write or append DataFrame to CSV.

    If append=True:
      - Appends rows
      - Writes header only if file does not exist

    If append=False:
      - Overwrites file
    """
    write_header = not (append and os.path.exists(filename))

    df.to_csv(
        filename,
        mode="a" if append else "w",
        header=write_header,
        index=False
    )





def dump_df_unique(
    df: pd.DataFrame,
    filename: str = "ED_data.csv",
    unique_cols: Optional[List[str]] = ["Lx", "Ly", "Nl", "beta", "J0", "J1", "J2", "h", "n_basis_vecs"],
):
    """
    Append DataFrame rows to CSV only if they are unique.

    Parameters
    ----------
    df : pd.DataFrame
        Data to append
    filename : str
        Target CSV file
    unique_cols : list[str] or None
        Columns defining uniqueness.
        If None, full-row uniqueness is used.
    """

    if not os.path.exists(filename):
        # File does not exist → write everything
        df.to_csv(filename, index=False)
        return

    # Load existing data
    existing = pd.read_csv(filename)

    if unique_cols is None:
        # Full-row uniqueness
        combined = pd.concat([existing, df], ignore_index=True)
        combined = combined.drop_duplicates()
    else:
        # Parameter-based uniqueness
        combined = pd.concat([existing, df], ignore_index=True)
        combined = combined.drop_duplicates(subset=unique_cols)

    # Overwrite file with deduplicated data
    combined.to_csv(filename, index=False)
