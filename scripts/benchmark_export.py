# Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
"""Export a completed local benchmark into the ignored Overleaf project folder."""
import argparse
from pathlib import Path
import sys
import json
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("identity", help="Benchmark job ID from the dashboard")
    parser.add_argument("--api", default="http://127.0.0.1:8000/api/benchmarks")
    args = parser.parse_args()
    with urllib.request.urlopen(f"{args.api}/{args.identity}", timeout=30) as response:
        job = json.load(response)
    if job["status"] != "complete":
        raise SystemExit("Complete the experiment before creating its publication bundle.")
    from app.services.benchmark import csv_result, figure_result
    import matplotlib as mpl
    mpl.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "font.family": "serif", "mathtext.fontset": "stix"})
    folder = ROOT / "docs/manuscript/benchmark" / args.identity
    for name in ("figures", "tables", "metadata"):
        (folder / name).mkdir(parents=True, exist_ok=True)
    (folder / "tables/observables.csv").write_text(csv_result(job), encoding="utf-8")
    (folder / "metadata/experiment.json").write_text(json.dumps(job, indent=2, allow_nan=False), encoding="utf-8")
    count = 0
    for axis in job["request"]["axes"]:
        for observable in job["observables"]:
            if observable.startswith("ct"):
                continue  # c(tau) is exported as a separate multi-time figure below.
            for format in ("pdf", "svg", "png"):
                (folder / f"figures/{axis}-{observable}.{format}").write_bytes(figure_result(job, axis, observable, format))
            count += 1
        for quantity in ("difference", "tau_int", "ess_per_second"):
            for observable in ("E", "ct0", "stripe2"):
                for format in ("pdf", "svg", "png"):
                    (folder / f"figures/{axis}-{observable}-{quantity}.{format}").write_bytes(figure_result(job, axis, observable, format, quantity))
                count += 1
    # Imaginary-time curves retain their physical time coordinate, separately
    # from Markov-chain ACFs and integrated correlation times.
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from app.services.benchmark import METHODS
    for point in job["points"]:
        figure = Figure(figsize=(6.4, 4.4), layout="constrained")
        figure.get_layout_engine().set(rect=(0, .065, 1, .935))
        FigureCanvasAgg(figure)
        ax = figure.subplots()
        ax.spines[["top", "right"]].set_visible(False)
        times = [0, .125, .25, .375, .5]
        ax.plot(times, [point["ed"]["values"][f"ct{i}"] for i in range(5)], "kD--", label="Complete ED (sampling SE = 0)", markersize=4)
        for method, color, marker in zip(METHODS, ("#0072B2", "#D55E00", "#009E73", "#CC79A7"), ("o", "s", "^", "v")):
            summary = point["methods"].get(method)
            if summary:
                for i, time in enumerate(times):
                    row = summary["observables"][f"ct{i}"]
                    ax.errorbar(time, row["mean"], yerr=row["error"], fmt=marker if row["status"] != "review" else "x", color=color, capsize=3)
                ax.plot([], [], marker, color=color, label=METHODS[method])
        ax.set(xlabel=r"$\tau/\beta$", ylabel=r"$c(\tau)$", title=f'{point["params"]["Lx"]} x {point["params"]["Ly"]}, '+", ".join(f"{k}={point['params'][k]:g}" for k in ("beta", "J1", "J2", "hx")))
        ax.legend(frameon=False, fontsize=7)
        figure.text(.02, .015, "QMC: 1 SE. Crosses require convergence review. This is imaginary time, not Markov-chain lag.", fontsize=6)
        for format in ("pdf", "svg", "png"):
            figure.savefig(folder / f'figures/correlation-{point["axis"]}-{point["x"]:g}.{format}', dpi=300)
        count += 1
    # A minimal standalone Overleaf entry point; no zip and no generated claims.
    lines = [r"\documentclass[aps,prb,twocolumn]{revtex4-2}", r"\usepackage{graphicx}", r"\begin{document}",
             r"\title{Finite-system ED and SSE-QMC validation}", r"\author{Vikas Vijigiri}", r"\maketitle",
             r"Complete ED references have zero sampling error. QMC bars denote one standard error estimated from correlated traces, blocked jackknife and independent chains. Cross markers indicate results requiring convergence review. This small-system comparison does not establish thermodynamic phase boundaries.",
             r"\section{Parameter sweeps}"]
    for axis in job["request"]["axes"]:
        for observable in ("E", "mx", "stripe2", "rho", "cv", "polar2"):
            lines += [r"\begin{figure}[ht]", rf"\includegraphics[width=\columnwidth]{{figures/{axis}-{observable}.pdf}}", rf"\caption{{Finite-system comparison: {axis.replace('_', ' ')} sweep, {observable.replace('_', ' ')}. Fixed parameters are recorded above each panel.}}", r"\end{figure}"]
    lines += [r"\end{document}"]
    (folder / "main.tex").write_text("\n".join(lines)+"\n", encoding="utf-8")
    print(f"Exported {count} figure panels in PDF, SVG and 300-dpi PNG to {folder}")


if __name__ == "__main__":
    main()
