"""Independent exact traces and quantum reference checks for the audit probes."""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import argparse
import itertools
import json
import subprocess
import numpy as np
from scipy.linalg import eigh
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh
import physics_audit as audit


def sparse_reference():
    energy, _, mz = audit.energies(4, 4, 1, .3, 0)
    rows = []
    for k in (20, 40, 64):
        w, v = eigsh(diags(energy), k=k, which="SA", v0=np.random.default_rng(10).normal(size=len(energy)), tol=1e-10)
        weights = np.exp(-10 * (w - w.min()))
        weights /= weights.sum()
        row = {"k": k, "ground_eigenpairs": int(np.sum(np.isclose(w, energy.min()))), "m2": float((mz * mz) @ (v * v) @ weights), "m4": float((mz ** 4) @ (v * v) @ weights), "E": float(w @ weights / 16), "max_residual": float(np.max(np.linalg.norm(energy[:, None] * v - v * w, axis=0)))}
        rows.append(row)
        print(json.dumps(row), flush=True)
    (audit.OUT / "sparse_ed.json").write_text(json.dumps(rows, indent=2))


def quantum_reference(intended=False):
    rows = []
    executable = audit.OUT / ("sse_probe_pair_bonds.exe" if intended else "sse_probe.exe")
    for j2, beta in ([(0, 1), (0, 10)] if intended else [(0, 1), (0, 10), (.025, 1), (.025, 10)]):
        energy, _, mz = audit.energies(4, 2, 1, .3, j2, intended)
        indices = np.arange(256)
        matrix = np.diag(energy)
        for site in range(8):
            matrix[indices, indices ^ (1 << site)] = -.3
        if not intended:
            plaquettes, neighbors = audit.ed_terms(4, 2)
            original = audit.ed.H2SQ_hamiltonian(8, 1, .3, j2, .3, plaquettes, neighbors).toarray()
            assert np.max(np.abs(matrix - original)) < 1e-12
        w, v = eigh(matrix)
        weights = np.exp(-beta * (w - w.min()))
        weights /= weights.sum()
        probability = (v * v) @ weights
        mx = sum(np.sum(v * v[indices ^ (1 << site), :], axis=0) for site in range(8)) @ weights / 8
        exact = {"E": float((w / 8) @ weights), "E2": float((w / 8) ** 2 @ weights), "E4": float((w / 8) ** 4 @ weights), "m2": float((mz * mz) @ probability), "m4": float((mz ** 4) @ probability), "Mx": float(mx)}
        row = {"lattice": [4, 2, 1], "beta": beta, "J2": j2, "hx": .3, "exact": exact, "runs": []}
        if not intended:
            reported = audit.ed.estimate_all_observables(4, 2, 1, 1, .3, j2, .3, w, v, beta).iloc[0]
            row["ED_reported_m2"] = float(reported["SMag_square"])
        for seed in (1, 2):
            name = f"{'intended' if intended else 'current'}_quantum_b{beta}_j2{j2}_seed{seed}"
            command = [str(executable), str(beta), "1", ".3", str(j2), ".3", "50000", "105000", str(seed), "hot", str(audit.OUT / (name + ".csv")), "4", "2"]
            result = json.loads(subprocess.check_output(command))
            result["seed"] = seed
            row["runs"].append(result)
        rows.append(row)
        print(json.dumps(row), flush=True)
        filename = "intended_quantum.json" if intended else "quantum_results.json"
        (audit.OUT / filename).write_text(json.dumps(rows, indent=2))


def local_balance():
    results = {}
    states = list(itertools.product((-1, 1), repeat=4))
    for name in ("J0", "J1"):
        weights = []
        transition = np.zeros((16, 16))
        for i, state in enumerate(states):
            a, b = state[0] * state[2], state[1] * state[3]
            weight = (2 if a * b == 1 else 0) if name == "J0" else (4 if a == b == 1 else 0)
            weights.append(weight)
            if not weight:
                continue
            pairs = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))] if name == "J0" else [((0, 2), (1, 3))]
            for pairing in pairs:
                for flips in itertools.product((0, 1), repeat=2):
                    target = list(state)
                    for pair, flip in zip(pairing, flips):
                        if flip:
                            for site in pair:
                                target[site] *= -1
                    transition[i, states.index(tuple(target))] += 1 / len(pairs) / 4
        flux = np.asarray(weights)[:, None] * transition
        assert np.max(np.abs(flux - flux.T)) < 1e-14
        assert not any(transition[i, j] > 0 and weights[j] == 0 for i in range(16) for j in range(16))
        results[name] = {"max_flux_asymmetry": float(np.max(np.abs(flux - flux.T)))}
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("sparse", "quantum", "intended-quantum", "balance"))
    mode = parser.parse_args().mode
    {"sparse": sparse_reference, "quantum": quantum_reference, "intended-quantum": lambda: quantum_reference(True), "balance": local_balance}[mode]()
