# Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
"""Finite-system reference for cyclic plaquettes and signed molecular dipoles."""
from functools import lru_cache
import numpy as np
from scipy.linalg import eigh

MODEL_VERSION = "cyclic-pair-j1-sublattice-dipoles-v3"
CX = np.array([1, -1, -1, 1]) / 4
CY = np.array([1, 1, -1, -1]) / 4
FRACTIONS = [0, .125, .25, .375, .5]


def geometry(lx, ly):
    if lx < 2 or ly < 2 or lx % 2 or ly % 2:
        raise ValueError("Periodic checkerboard dimensions must be even and at least two.")
    positions = [(x, y) for y in range(ly) for x in range(lx) if (x+y) % 2 == 0]
    index = {position: i for i, position in enumerate(positions)}
    plaquettes = [[x+lx*y, x+lx*((y+1) % ly), (x+1) % lx+lx*((y+1) % ly), (x+1) % lx+lx*y] for x, y in positions]
    # Periodic bonds are a multiset: coincident directions at L=2 remain distinct.
    pairs = [(i, index[((x+dx) % lx, (y+1) % ly)]) for i, (x, y) in enumerate(positions) for dx in (1, -1)]
    signs = np.array([1 if x % 2 == 0 else -1 for x, _ in positions])
    return np.asarray(plaquettes), pairs, signs


@lru_cache(maxsize=16)
def basis_data(lx, ly):
    n = lx*ly
    if n > 16:
        raise ValueError("Complete diagonal references are limited to 16 spins.")
    indices = np.arange(1 << n)
    spins = 2*((indices[:, None] >> np.arange(n)) & 1)-1
    plaquettes, pairs, signs = geometry(lx, ly)
    local = spins[:, plaquettes]
    px, py = local @ CX, local @ CY
    oriented_x, oriented_y = px*signs, py*signs
    parity = np.prod(local, axis=2).mean(axis=1)
    ice = (local.sum(axis=2) == 0).mean(axis=1)
    x, y = np.arange(n) % lx, np.arange(n) // lx
    mz = spins @ ((-1.)**(x+y)) / (2*n)
    stripe2 = (spins @ ((-1.)**x) / (2*n))**2 + (spins @ ((-1.)**y) / (2*n))**2
    polar2 = oriented_x.mean(axis=1)**2 + oriented_y.mean(axis=1)**2
    observables = {"mz": mz, "mz2": mz**2, "mz4": mz**4, "plaquette": parity, "ice": ice, "rho": (4*ice-1)/3, "stripe2": stripe2, "stripe4": stripe2**2, "polar2": polar2, "polar4": polar2**2}
    return indices, local, px, py, pairs, observables, oriented_x.mean(axis=1), oriented_y.mean(axis=1)


@lru_cache(maxsize=24)
def spectrum(lx, ly, j0, j1, j2, hx):
    indices, local, px, py, pairs, observables, polx, poly = basis_data(lx, ly)
    diagonal = -j0*np.prod(local, axis=2).sum(axis=1)
    diagonal = diagonal + j1*(local[:, :, 0]*local[:, :, 2] + local[:, :, 1]*local[:, :, 3]).sum(axis=1)
    # P_B=-tilde(P_B), so -J2 P_A.P_B = +J2 tilde(P_A).tilde(P_B).
    for a, b in pairs:
        diagonal = diagonal + j2*(px[:, a]*px[:, b] + py[:, a]*py[:, b])
    if hx == 0:
        return diagonal, None, diagonal, 0.0
    if lx*ly > 10:
        raise ValueError("Complete quantum ED supports at most 10 spins; use 4x2, or h=0 for 4x4.")
    matrix = np.diag(diagonal)
    for site in range(lx*ly):
        matrix[indices, indices ^ (1 << site)] = -hx
    w, v = eigh(matrix, driver="evd", check_finite=False)
    residual = float(np.max(np.linalg.norm(matrix @ v-v*w, axis=0)))
    return w, v, diagonal, residual


def exact_reference(params):
    lx, ly = int(params["Lx"]), int(params["Ly"])
    n = lx*ly
    beta, hx = float(params["beta"]), float(params.get("hx", 0))
    indices, local, px, py, pairs, diagonals, polx, poly = basis_data(lx, ly)
    w, v, diagonal, residual = spectrum(lx, ly, float(params["J0"]), float(params["J1"]), float(params["J2"]), hx)
    shifted = w-w.min()
    boltzmann = np.exp(-beta*shifted)
    weights = boltzmann / boltzmann.sum()
    probabilities = weights if v is None else (v*v) @ weights
    values = {key: float(value @ probabilities) for key, value in diagonals.items()}
    values.update(E=float((w/n) @ weights), E2=float((w/n)**2 @ weights), E4=float((w/n)**4 @ weights))
    values["cv"] = max(0., beta*beta*n*(values["E2"]-values["E"]**2))
    values["binder"] = .5*(1-values["stripe4"]/(3*values["stripe2"]**2)) if values["stripe2"] > 0 else None
    ct = []
    if v is None:
        values["mx"] = 0.0
        values["chi_stag"] = beta*n*(values["mz2"]-values["mz"]**2)
        values["chi_polar"] = beta*len(pairs)/2*((polx*polx+poly*poly) @ weights - (polx @ weights)**2-(poly @ weights)**2)
        delta = diagonal[indices[:, None] ^ (1 << np.arange(n))]-diagonal[:, None]
        other_weights = weights[indices[:, None] ^ (1 << np.arange(n))]
        kernel = np.divide(weights[:, None]-other_weights, delta, out=np.broadcast_to(beta*weights[:, None], delta.shape).copy(), where=np.abs(delta) > 1e-12)
        values["chi_x"] = float(kernel.sum()/n)
        ct = [float((px*px).mean(axis=1) @ weights)]*len(FRACTIONS)
    else:
        difference = w[None, :]-w[:, None]
        kubo = np.divide(weights[:, None]-weights[None, :], difference, out=np.broadcast_to(beta*weights[:, None], difference.shape).copy(), where=np.abs(difference) > 1e-12)
        def response(operator, normalization):
            rotated = v.T @ (operator[:, None]*v)
            return float(normalization*((kubo*rotated**2).sum()-beta*(operator @ probabilities)**2))
        values["chi_stag"] = response(diagonals["mz"], n)
        values["chi_polar"] = response(polx, local.shape[1]) + response(poly, local.shape[1])
        xv = sum(v[indices ^ (1 << site)] for site in range(n))/n
        rotated_x = v.T @ xv
        values["mx"] = float(np.diag(rotated_x) @ weights)
        values["chi_x"] = float(n*((kubo*rotated_x**2).sum()-beta*values["mx"]**2))
        rotated_dipoles = [v.T @ (px[:, p, None]*v) for p in range(local.shape[1])]
        amplitudes = sum(rotated**2 for rotated in rotated_dipoles)/local.shape[1]
        for fraction in FRACTIONS:
            thermal = np.exp(-beta*((1-fraction)*shifted[:, None]+fraction*shifted[None, :])) / boltzmann.sum()
            ct.append(float((thermal*amplitudes).sum()))
    values.update({f"ct{i}": value for i, value in enumerate(ct)})
    return {"values": values, "statistical_errors": {key: 0.0 for key in values}, "numerical_residual": residual, "reference": "complete classical trace" if v is None else "complete dense eigenspectrum", "dimension": len(w), "ground_energy_per_spin": float(w.min()/n), "model_version": MODEL_VERSION}
