# Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
"""Export a completed local benchmark into the ignored Overleaf project folder."""
import argparse
from pathlib import Path
import sys
import json
import re
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("identity", help="Benchmark job ID from the dashboard")
    parser.add_argument("--name", help="Optional publication folder label (letters, digits, underscores and hyphens)")
    parser.add_argument("--api", default="http://127.0.0.1:8000/api/benchmarks")
    args = parser.parse_args()
    if args.name and not re.fullmatch(r"[A-Za-z0-9_-]+", args.name):
        raise SystemExit("Use letters, digits, underscores and hyphens for the folder label.")
    with urllib.request.urlopen(f"{args.api}/{args.identity}", timeout=30) as response:
        job = json.load(response)
    if job["status"] != "complete":
        raise SystemExit("Complete the experiment before creating its publication bundle.")
    from app.services.benchmark import csv_result, figure_result
    import matplotlib as mpl
    mpl.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "font.family": "serif", "mathtext.fontset": "stix"})
    folder = ROOT / "docs/manuscript/benchmark" / (args.name or args.identity)
    existing=folder/"metadata/experiment.json"
    if existing.exists():
        previous=json.loads(existing.read_text(encoding="utf-8")).get("id")
        if previous!=job["id"]:
            label=re.sub(r"[^A-Za-z0-9_-]","",str(previous))[:32] or "previous"
            archived=ROOT/"archive/benchmark-exports"/f"{folder.name}-{label}"
            suffix=1
            while archived.exists():
                archived=archived.with_name(f"{folder.name}-{label}-{suffix}");suffix+=1
            if not folder.resolve().is_relative_to(ROOT) or not archived.resolve().is_relative_to(ROOT):
                raise SystemExit("Publication archive paths must remain inside this workspace.")
            archived.parent.mkdir(parents=True,exist_ok=True)
            folder.rename(archived)
    for name in ("figures", "tables", "metadata"):
        (folder / name).mkdir(parents=True, exist_ok=True)
    (folder / "tables/observables.csv").write_text(csv_result(job), encoding="utf-8")
    (folder / "metadata/experiment.json").write_text(json.dumps(job, indent=2, allow_nan=False), encoding="utf-8")
    import csv
    with (folder/"tables/precision.csv").open("w",newline="",encoding="utf-8") as output:
        writer=csv.writer(output);writer.writerow(["axis","x","schedule","observable","runs","measurement_rows","mean","SE","measured_SE_ratio","expected_SE_ratio"])
        for point in job["points"]:
            for method,summary in point["methods"].items():
                for row in summary.get("error_scaling",[]):
                    for key,observation in row["observables"].items():
                        writer.writerow([point["axis"],point["x"],method,key,row["chains"],row["samples"],observation["mean"],observation["error"],observation.get("se_ratio"),row["expected_se_ratio"]])
    with (folder/"tables/estimator-gains.csv").open("w",newline="",encoding="utf-8") as output:
        writer=csv.writer(output);writer.writerow(["axis","x","schedule","observable","boundary_mean","boundary_SE","improved_mean","improved_SE","SE_gain"])
        for point in job["points"]:
            for method,summary in point["methods"].items():
                for key,row in summary.get("estimator_comparison",{}).items():
                    writer.writerow([point["axis"],point["x"],method,key,*[row[k] for k in ["boundary_mean","boundary_SE","improved_mean","improved_SE","SE_gain"]]])
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
        for method, color, marker in zip(METHODS, ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00"), ("o", "s", "^", "v", "P")):
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
    precision_names=[]
    seen=set()
    for point in job["points"]:
        if point["point_id"] in seen:continue
        seen.add(point["point_id"])
        if not any(summary.get("error_scaling") for summary in point["methods"].values()):continue
        for observable in ("E","stripe2","ct0"):
            figure=Figure(figsize=(6.4,6.2),layout="constrained")
            figure.get_layout_engine().set(rect=(0,.05,1,.95));FigureCanvasAgg(figure)
            axes=figure.subplots(2,1)
            for ax in axes:
                ax.spines[["top","right"]].set_visible(False);ax.grid(alpha=.15)
                ax.set_xlabel("Independent runs per schedule")
            axes[0].axhline(point["ed"]["values"][observable],color="black",linestyle="--",label="Complete ED")
            expected_drawn=False
            for method,color,marker in zip(METHODS,("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00"),("o","s","^","v","P")):
                history=point["methods"].get(method,{}).get("error_scaling",[])
                valid=[row for row in history if row["observables"][observable].get("error") is not None]
                if not valid:continue
                x=[row["chains"] for row in valid]
                axes[0].errorbar(x,[r["observables"][observable]["mean"] for r in valid],yerr=[r["observables"][observable]["error"] for r in valid],fmt=marker+"-",color=color,capsize=3,markersize=4,label=METHODS[method])
                axes[1].plot(x,[r["observables"][observable].get("se_ratio") for r in valid],marker+"-",color=color,markersize=4)
                axes[1].plot(x,[r["expected_se_ratio"] for r in valid],"--",color=color,alpha=.65,label="Expected square-root scaling" if not expected_drawn else None)
                expected_drawn=True
            axes[0].set_ylabel(job["observables"][observable][1]);axes[1].set_ylabel("SE / first-batch SE")
            axes[0].legend(frameon=False,fontsize=6);axes[1].legend(frameon=False,fontsize=7)
            figure.suptitle(f'{point["params"]["Lx"]} x {point["params"]["Ly"]} periodic | '+", ".join(f"{k}={point['params'][k]:g}" for k in ("beta","J1","J2","hx")),fontsize=9)
            figure.text(.02,.012,"QMC bars: 1 SE. Observed errors may fluctuate; shrinking error alone does not establish correct sampling.",fontsize=6)
            name=f"precision-{point['axis']}-{point['x']:g}-{observable}"
            for format in ("pdf","svg","png"):figure.savefig(folder/f"figures/{name}.{format}",dpi=300)
            precision_names.append(name);count+=1
    # A minimal standalone Overleaf entry point; no zip and no generated claims.
    lines = [r"\documentclass[aps,prb,twocolumn]{revtex4-2}", r"\usepackage{graphicx}", r"\begin{document}",
             r"\title{Finite-system ED and SSE-QMC validation}", r"\author{Vikas Vijigiri}", r"\maketitle",
             r"Complete ED references have zero sampling error. QMC bars denote one standard error estimated from correlated traces, blocked jackknife and independent chains. Cross markers indicate results requiring convergence review. This small-system comparison does not establish thermodynamic phase boundaries.",
             r"\section{Parameter sweeps}"]
    for axis in job["request"]["axes"]:
        for observable in ("E", "mx", "stripe2", "rho", "cv", "polar2"):
            lines += [r"\begin{figure}[ht]", rf"\includegraphics[width=\columnwidth]{{figures/{axis}-{observable}.pdf}}", rf"\caption{{Finite-system comparison: {axis.replace('_', ' ')} sweep, {observable.replace('_', ' ')}. Fixed parameters are recorded above each panel.}}", r"\end{figure}"]
    if precision_names:
        lines += [r"\section{Statistical precision}",r"Independent trajectories are pooled within each update schedule. Error estimates account for correlations and independent-run dispersion. The dashed curves show square-root sample-size scaling; actual estimates need not decrease at every intermediate batch. Improved measurements average propagated SSE states and use exact Kubo sums; see A. W. Sandvik, arXiv:1101.3281 and arXiv:1909.10591."]
        for name in precision_names[:3]:
            lines += [r"\begin{figure}[ht]",rf"\includegraphics[width=\columnwidth]{{figures/{name}.pdf}}",r"\caption{Cumulative estimates and measured standard-error scaling.}",r"\end{figure}"]
    lines += [r"\end{document}"]
    (folder / "main.tex").write_text("\n".join(lines)+"\n", encoding="utf-8")
    print(f"Exported {count} figure panels in PDF, SVG and 300-dpi PNG to {folder}")


if __name__ == "__main__":
    main()
