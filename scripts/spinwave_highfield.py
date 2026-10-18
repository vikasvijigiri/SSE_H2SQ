# Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
"""High-field harmonic expansion checked against exact small-stack spectra.

The expansion uses the Pauli Hamiltonian directly.  Its output is a local
cross-check of the x-polarized high-field branch, not a replacement for the
low-field spin-wave theory missing from the historical source tree.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))
from app.services.physics_model import geometry  # noqa: E402

OUTPUT = ROOT / ".cache" / "spinwave-highfield"
PAPER_FIGURES = ROOT / "manuscript_latest" / "review" / "overleaf" / "figures"
PAPER_TABLES = ROOT / "manuscript_latest" / "review" / "overleaf" / "tables"


def pair_couplings(lx: int, ly: int, lz: int, j1: float, j2: float, j3: float) -> dict[tuple[int, int], float]:
    """Expand all two-spin terms into unordered pairs; discard constants."""
    plaquettes, neighbors, _ = geometry(lx, ly)
    area = lx * ly
    edges: dict[tuple[int, int], float] = {}

    def add(i: int, j: int, value: float) -> None:
        if i == j:
            return  # sigma_i^2 is a constant and has no harmonic curvature.
        key = (min(i, j), max(i, j))
        edges[key] = edges.get(key, 0.0) + float(value)

    for z in range(lz):
        off = z * area
        for p in plaquettes:
            add(off + int(p[0]), off + int(p[2]), j1)
            add(off + int(p[1]), off + int(p[3]), j1)
        # The manuscript has -J2 P_A.P_B.  Since P_B=-Ptilde_B,
        # this is +J2 Ptilde_A.Ptilde_B, as in the ED/SSE implementation.
        cx = (1.0, -1.0, -1.0, 1.0)
        cy = (1.0, 1.0, -1.0, -1.0)
        for a, b in neighbors:
            for u in range(4):
                for v in range(4):
                    value = j2 * (cx[u] * cx[v] + cy[u] * cy[v]) / 16.0
                    add(off + int(plaquettes[a, u]), off + int(plaquettes[b, v]), value)
    for z in range(lz):
        zn = (z + 1) % lz
        for site in range(area):
            add(z * area + site, zn * area + site, j3)
    return {edge: value for edge, value in edges.items() if abs(value) > 1e-14}


def diagonal_energy(lx: int, ly: int, lz: int, j0: float, j1: float, j2: float, j3: float) -> np.ndarray:
    """Independent diagonal evaluation in the sigma-z basis."""
    area = lx * ly
    n = area * lz
    plaquettes, neighbors, _ = geometry(lx, ly)
    states = 2 * ((np.arange(1 << n)[:, None] >> np.arange(n)) & 1) - 1
    total = np.zeros(len(states), dtype=float)
    for z in range(lz):
        layer = states[:, z * area:(z + 1) * area]
        local = layer[:, plaquettes]
        total -= j0 * np.prod(local, axis=2).sum(axis=1)
        total += j1 * (local[:, :, 0] * local[:, :, 2] + local[:, :, 1] * local[:, :, 3]).sum(axis=1)
        px = local @ (np.array([1., -1., -1., 1.]) / 4.)
        py = local @ (np.array([1., 1., -1., -1.]) / 4.)
        for a, b in neighbors:
            total += j2 * (px[:, a] * px[:, b] + py[:, a] * py[:, b])
    for z in range(lz):
        total += j3 * np.sum(states[:, z * area:(z + 1) * area] *
                             states[:, ((z + 1) % lz) * area:((z + 1) % lz + 1) * area], axis=1)
    return total


def expanded_pair_energy(lx: int, ly: int, lz: int, j0: float, j1: float, j2: float, j3: float) -> np.ndarray:
    """Pair-graph expansion plus the four-spin plaquette term."""
    area = lx * ly
    n = area * lz
    plaquettes, _, _ = geometry(lx, ly)
    states = 2 * ((np.arange(1 << n)[:, None] >> np.arange(n)) & 1) - 1
    result = np.zeros(len(states), dtype=float)
    for z in range(lz):
        result -= j0 * np.prod(states[:, z * area:(z + 1) * area][:, plaquettes], axis=2).sum(axis=1)
    for (i, j), value in pair_couplings(lx, ly, lz, j1, j2, j3).items():
        result += value * states[:, i] * states[:, j]
    return result


def harmonic_modes(lx: int, ly: int, lz: int, j1: float, j2: float, j3: float, field: float, spin: float = 0.5) -> np.ndarray:
    """Positive frequencies of the x-polarized quadratic HP Hamiltonian."""
    n = lx * ly * lz
    adjacency = np.zeros((n, n), dtype=float)
    for (i, j), value in pair_couplings(lx, ly, lz, j1, j2, j3).items():
        adjacency[i, j] += value
        adjacency[j, i] += value
    lam = np.linalg.eigvalsh(adjacency)
    radicand = field * (field + 2.0 * spin * lam)
    if np.min(radicand) < -1e-12:
        return np.full(n, np.nan)
    return 2.0 * np.sqrt(np.maximum(radicand, 0.0))


def exact_odd_gap(lx: int, ly: int, lz: int, j0: float, j1: float, j2: float, j3: float, field: float) -> float:
    """Lowest odd-global-x-parity gap from exact finite-system diagonalization."""
    n = lx * ly * lz
    if n > 10:
        raise ValueError("The dense parity-resolved reference is limited to 10 spins.")
    diag = diagonal_energy(lx, ly, lz, j0, j1, j2, j3)
    dim = 1 << n
    matrix = np.diag(diag)
    indices = np.arange(dim)
    for site in range(n):
        matrix[indices, indices ^ (1 << site)] = -field
    representatives = np.arange(dim // 2)
    complement = (dim - 1) - representatives
    gaps = {}
    for parity in (+1, -1):
        basis = np.zeros((dim, dim // 2), dtype=float)
        cols = np.arange(dim // 2)
        basis[representatives, cols] = 1.0 / np.sqrt(2.0)
        basis[complement, cols] = parity / np.sqrt(2.0)
        block = basis.T @ matrix @ basis
        gaps[parity] = float(eigh(block, eigvals_only=True, subset_by_index=[0, 0], check_finite=False)[0])
    return gaps[-1] - min(gaps.values())


def bloch_bond_shells(j1: float, j2: float, j3: float,
                      lx: int = 8, ly: int = 8, lz: int = 4) -> list[list[tuple[int, int, int, int, float]]]:
    """Extract directed pair bonds from one A/B spin basis in the checkerboard cell.

    The checkerboard is periodic under a1=(1,1), a2=(1,-1).  A and B are the
    even- and odd-parity physical spin sites, represented by (0,0) and (1,0).
    A larger periodic cell avoids the coincident-bond special cases at L=2.
    Entries are (other_sublattice, dx, dy, dz, J_ij), using physical-site
    displacement so the Fourier phase is independent of basis gauge.
    """
    area = lx * ly
    edges = pair_couplings(lx, ly, lz, j1, j2, j3)
    shells: list[list[tuple[int, int, int, int, float]]] = [[], []]
    roots = ((0, 0), (1, 0))

    def minimum_image(delta: int, length: int) -> int:
        return (delta + length // 2) % length - length // 2

    for alpha, (x, y) in enumerate(roots):
        root = x + lx * y + area  # take z=1, away from the vertical boundary
        for (i, j), coupling in edges.items():
            if i == root:
                other = j
            elif j == root:
                other = i
            else:
                continue
            z, rem = divmod(other, area)
            oy, ox = divmod(rem, lx)
            dx = minimum_image(ox - x, lx)
            dy = minimum_image(oy - y, ly)
            dz = minimum_image(z - 1, lz)
            beta = (ox + oy) % 2
            shells[alpha].append((beta, dx, dy, dz, coupling))
    return shells


def bloch_adjacency(q1: float, q2: float, kz: float,
                    shells: list[list[tuple[int, int, int, int, float]]]) -> np.ndarray:
    """2x2 Hermitian pair-coupling matrix at reduced momentum (q1,q2,kz).

    q1=k dot a1 and q2=k dot a2 for the rotated checkerboard primitive cell.
    The corresponding physical in-plane wave vector is
    kx=(q1+q2)/2, ky=(q1-q2)/2 in inverse square-lattice spacings.
    """
    kx = 0.5 * (q1 + q2)
    ky = 0.5 * (q1 - q2)
    matrix = np.zeros((2, 2), dtype=np.complex128)
    for alpha, bonds in enumerate(shells):
        for beta, dx, dy, dz, coupling in bonds:
            phase = np.exp(1j * (kx * dx + ky * dy + kz * dz))
            matrix[alpha, beta] += coupling * phase
    if not np.allclose(matrix, matrix.conj().T, atol=1e-11, rtol=0.):
        raise ArithmeticError("Bloch pair-coupling matrix is not Hermitian.")
    return matrix


def spinwave_bands(qpoints: np.ndarray, j1: float, j2: float, j3: float,
                   field: float, spin: float = 0.5) -> np.ndarray:
    """Harmonic HP bands at rows (q1,q2,kz), in the same units as the couplings."""
    shells = bloch_bond_shells(j1, j2, j3)
    eigenvalues = np.asarray([
        np.linalg.eigvalsh(bloch_adjacency(q1, q2, kz, shells))
        for q1, q2, kz in np.asarray(qpoints, dtype=float)
    ])
    radicand = field * (field + 2.0 * spin * eigenvalues)
    if np.any(radicand < -1e-12):
        raise ValueError("The fully x-polarized state is harmonically unstable on this path.")
    return 2.0 * np.sqrt(np.maximum(radicand, 0.0))


def high_symmetry_path(points_per_segment: int = 80) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    """Simple-tetragonal path in q=(k.a1,k.a2,kz), returned in radians."""
    if points_per_segment < 2:
        raise ValueError("Use at least two points per high-symmetry segment.")
    vertices = np.pi * np.asarray([
        (0., 0., 0.),  # Gamma
        (1., 0., 0.),  # X
        (1., 1., 0.),  # M
        (0., 0., 0.),  # Gamma
        (0., 0., 1.),  # Z
        (1., 0., 1.),  # R
        (1., 1., 1.),  # A
        (0., 0., 1.),  # Z
    ])
    labels = [r"$\Gamma$", "X", "M", r"$\Gamma$", "Z", "R", "A", "Z"]
    qpoints = []
    ticks = [0.]
    for start, stop in zip(vertices[:-1], vertices[1:]):
        qpoints.extend(start + t * (stop - start)
                       for t in np.linspace(0., 1., points_per_segment, endpoint=False))
        ticks.append(ticks[-1] + np.linalg.norm(_reduced_to_physical(stop) - _reduced_to_physical(start)))
    qpoints.append(vertices[-1])
    qpoints = np.asarray(qpoints)
    physical = np.asarray([_reduced_to_physical(q) for q in qpoints])
    distance = np.concatenate(([0.], np.cumsum(np.linalg.norm(np.diff(physical, axis=0), axis=1))))
    return qpoints, distance, np.asarray(ticks), labels


def _reduced_to_physical(qpoint: np.ndarray) -> np.ndarray:
    q1, q2, kz = qpoint
    return np.asarray([0.5 * (q1 + q2), 0.5 * (q1 - q2), kz])


def run(lx: int = 2, ly: int = 2, lz: int = 2, j0: float = 1.0, j1: float = .3,
        j2: float = .025, j3_values: tuple[float, ...] = (-.3, 0., .3),
        fields: tuple[float, ...] = (2., 4., 8., 16.)) -> list[dict]:
    """Run the Pauli-factor, pair-expansion, and ED/HP gap cross-checks."""
    rows = []
    for j3 in j3_values:
        direct = diagonal_energy(lx, ly, lz, j0, j1, j2, j3)
        expanded = expanded_pair_energy(lx, ly, lz, j0, j1, j2, j3)
        offset = float(np.mean(direct - expanded))
        expansion_residual = float(np.max(np.abs(direct - expanded - offset)))
        for field in fields:
            modes = harmonic_modes(lx, ly, lz, j1, j2, j3, field)
            ed = exact_odd_gap(lx, ly, lz, j0, j1, j2, j3, field)
            hp = float(np.min(modes))
            rows.append({"J3": j3, "field": field, "ED_odd_gap": ed,
                         "HP_min_gap": hp, "abs_difference": abs(ed - hp),
                         "relative_difference": abs(ed - hp) / ed if ed else None,
                         "pair_expansion_constant": offset,
                         "pair_expansion_max_residual": expansion_residual,
                         "HP_stable": bool(np.all(np.isfinite(modes)))})
    return rows


def save(rows: list[dict], outdir: Path = OUTPUT) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    with (outdir / "highfield-gap-check.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (outdir / "highfield-gap-check.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")
    PAPER_TABLES.mkdir(parents=True, exist_ok=True)
    table = [
        r"\begin{table*}[t]", r"\caption{Lowest odd-global-$x$-parity excitation gap on the periodic $2\times2\times2$ stack. ED is exact within the eight-spin Hilbert space; the harmonic column is the minimum high-field Holstein--Primakoff frequency. The relative difference is $|\Delta_{\rm HP}-\Delta_{\rm ED}|/\Delta_{\rm ED}$.}",
        r"\label{tab:highfield-gap}", r"\begin{ruledtabular}",
        r"\begin{tabular}{rrrrr}", r"$J_3/J_0$ & $K/J_0$ & $\Delta_{\rm ED}/J_0$ & $\Delta_{\rm HP}/J_0$ & Relative difference \\", r"\hline"
    ]
    for row in sorted(rows, key=lambda item: (item["J3"], item["field"])):
        table.append(f'{row["J3"]:+.1f} & {row["field"]:.1f} & {row["ED_odd_gap"]:.6f} & {row["HP_min_gap"]:.6f} & {row["relative_difference"]:.4f} \\\\')
    table.extend([r"\end{tabular}", r"\end{ruledtabular}", r"\end{table*}"])
    (PAPER_TABLES / "highfield-gap-check.tex").write_text("\n".join(table) + "\n", encoding="utf-8")

    plt.rcParams.update({"font.family": "serif", "font.serif": ["STIX Two Text", "STIXGeneral", "DejaVu Serif"],
                         "mathtext.fontset": "stix", "font.size": 8.2, "axes.labelsize": 8.2,
                         "axes.titlesize": 8.5, "legend.fontsize": 7.3, "xtick.labelsize": 7.5,
                         "ytick.labelsize": 7.5, "axes.linewidth": .7, "lines.linewidth": 1.25,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    # The top panels show the excitation scale; the lower panels expose the
    # high-field convergence that is invisible when ED and HP nearly overlap.
    fig, axes = plt.subplots(2, 2, figsize=(7.15, 4.55), sharex="col",
                             gridspec_kw={"height_ratios": [2.1, 1.0]})
    colors = {"ED": "#263B53", "HP": "#D55E00"}
    for col, j3 in enumerate((0.0, 0.3)):
        selected = sorted((r for r in rows if r["J3"] == j3), key=lambda row: row["field"])
        x = np.array([r["field"] for r in selected])
        ax = axes[0, col]
        ax.plot(x, [r["ED_odd_gap"] for r in selected], marker="s", ms=4.0,
                color=colors["ED"], label="Exact diagonalization")
        ax.plot(x, [r["HP_min_gap"] for r in selected], marker="o", ms=4.1,
                mfc="white", mew=1.0, linestyle="--", color=colors["HP"],
                label="Harmonic HP")
        panel = "(a)" if col == 0 else "(b)"
        title = (r"$J_3/J_0=0$" if j3 == 0 else r"$|J_3|/J_0=0.3$ (both signs coincide)")
        ax.set_title(f"{panel}  {title}", loc="left", pad=5)
        ax.set_xscale("log", base=2)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{value:g}" for value in x])
        ax.tick_params(direction="in", top=True, right=True, length=3.0, width=.7)
        ax.margins(x=.08, y=.12)
        for spine in ax.spines.values():
            spine.set_linewidth(.65)
        error_ax = axes[1, col]
        rel = 100.0 * np.array([r["relative_difference"] for r in selected])
        error_ax.plot(x, rel, marker="D", ms=3.8, color="#7A5195", lw=1.1)
        error_ax.set_xscale("log", base=2)
        error_ax.set_yscale("log")
        error_ax.set_xticks(x)
        error_ax.set_xticklabels([f"{value:g}" for value in x])
        error_ax.set_ylim(.1, 16.)
        error_ax.set_yticks([.1, 1., 10.], ["0.1", "1", "10"])
        for ref in (0.1, 1., 10.):
            error_ax.axhline(ref, color="#98A2B3", lw=.55, ls=(0, (2, 2)), zorder=0)
        error_ax.tick_params(direction="in", top=True, right=True, length=3.0, width=.7)
        error_ax.set_xlabel(r"Transverse field $K/J_0$")
        error_ax.margins(x=.08)
        error_ax.text(.98, .90, "(c)" if col == 0 else "(d)", transform=error_ax.transAxes,
                      ha="right", va="top", fontweight="bold", fontsize=8.2)
        for spine in error_ax.spines.values():
            spine.set_linewidth(.65)
    axes[0, 0].set_ylabel(r"Lowest odd-parity gap, $\Delta/J_0$")
    axes[1, 0].set_ylabel(r"$|\Delta_{\rm HP}-\Delta_{\rm ED}|/\Delta_{\rm ED}$ (%)")
    fig.legend(*axes[0, 0].get_legend_handles_labels(), frameon=False, ncol=2,
               loc="upper center", bbox_to_anchor=(.52, .995), handlelength=2.2,
               columnspacing=2.0)
    fig.subplots_adjust(left=.105, right=.985, bottom=.13, top=.86, hspace=.12, wspace=.22)
    for suffix in ("pdf", "svg", "png"):
        fig.savefig(outdir / f"highfield-gap-check.{suffix}", dpi=600 if suffix == "png" else None)
    PAPER_FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(PAPER_FIGURES / "highfield-gap-check.pdf")
    plt.close(fig)


def save_spectrum_figure(outdir: Path, j0: float, j1: float, j2: float,
                         j3_values: tuple[float, ...], field: float) -> None:
    """Save the two-sublattice high-symmetry spectrum for representative J3 values."""
    if j0 <= 0.:
        raise ValueError("J0 must be positive to plot dimensionless frequencies.")
    qpoints, distance, ticks, labels = high_symmetry_path()
    j3_values = tuple(j3_values)
    if not j3_values:
        raise ValueError("At least one J3 value is required for the spectrum plot.")
    plt.rcParams.update({"font.family": "serif", "font.serif": ["STIX Two Text", "STIXGeneral", "DejaVu Serif"],
                         "mathtext.fontset": "stix", "font.size": 8.3, "axes.labelsize": 8.5,
                         "axes.titlesize": 8.5, "legend.fontsize": 7.4, "xtick.labelsize": 7.8,
                         "ytick.labelsize": 7.8, "axes.linewidth": .7, "lines.linewidth": 1.45,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axs = plt.subplots(1, len(j3_values), figsize=(7.15, 3.15), squeeze=False,
                            sharey=True)
    colors = ("#0072B2", "#D55E00")
    for panel, (ax, j3) in enumerate(zip(axs[0], j3_values)):
        bands = spinwave_bands(qpoints, j1, j2, j3, field) / j0
        ax.plot(distance, bands[:, 0], color=colors[0], label="Lower branch")
        ax.plot(distance, bands[:, 1], color=colors[1], label="Upper branch")
        ax.set_title(fr"({chr(ord('a') + panel)})  $J_3/J_0={j3 / j0:g}$", loc="left", pad=6)
        ax.set_xticks(ticks, labels)
        for tick in ticks[1:-1]:
            ax.axvline(tick, color="#AAB3BE", lw=.55, ls=(0, (2, 2)), zorder=0)
        ax.tick_params(direction="in", top=True, right=True, length=3., width=.7)
        ax.margins(x=0.)
        for spine in ax.spines.values():
            spine.set_linewidth(.65)
        ax.set_xlabel("Wave vector along high-symmetry path")
    axs[0, 0].set_ylabel(r"Harmonic frequency, $\omega/J_0$")
    fig.suptitle(fr"High-field spin-wave spectrum: $K/J_0={field / j0:g}$, $J_1/J_0={j1 / j0:g}$, $J_2/J_0={j2 / j0:g}$",
                 y=.985, fontsize=8.5)
    handles, legend_labels = axs[0, 0].get_legend_handles_labels()
    fig.legend(handles, legend_labels, frameon=False, ncol=2, loc="upper center",
               bbox_to_anchor=(.5, .91), handlelength=2.0, columnspacing=2.0)
    fig.subplots_adjust(left=.095, right=.99, bottom=.22, top=.79, wspace=.13)
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / "highfield-spectrum.pdf")
    fig.savefig(outdir / "highfield-spectrum.png", dpi=600)
    PAPER_FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(PAPER_FIGURES / "highfield-spectrum.pdf")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--lx", type=int, default=2)
    parser.add_argument("--ly", type=int, default=2)
    parser.add_argument("--lz", type=int, default=2)
    parser.add_argument("--j0", type=float, default=1.)
    parser.add_argument("--j1", type=float, default=.3)
    parser.add_argument("--j2", type=float, default=.025)
    parser.add_argument("--j3", type=float, nargs="+", default=[-.3, 0., .3])
    parser.add_argument("--fields", type=float, nargs="+", default=[2., 4., 8., 16.])
    parser.add_argument("--spectrum-field", type=float, default=None,
                        help="field for the high-symmetry bands (defaults to 4 J0)")
    args = parser.parse_args()
    rows = run(args.lx, args.ly, args.lz, args.j0, args.j1, args.j2,
               tuple(args.j3), tuple(args.fields))
    save(rows, args.output)
    spectrum_field = 4. * args.j0 if args.spectrum_field is None else args.spectrum_field
    spectrum_j3 = tuple(sorted({0., max((abs(value) for value in args.j3), default=0.)}))
    save_spectrum_figure(args.output, args.j0, args.j1, args.j2, spectrum_j3, spectrum_field)
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
