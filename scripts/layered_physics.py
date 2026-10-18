# Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
"""Finite-stack ED/SSE benchmark for layer parity and local ice diagnostics.

The default study uses the periodic 2x2x2 checkerboard stack, which has only
eight pseudospins and therefore admits a complete dense thermal ED reference.
It checks the exact uniform-J3 sign mapping and resolves qz=0 and qz=pi dipole
structure factors. Outputs stay in the ignored local .cache tree; only source
and plotting code are intended for version control.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))
from app.services.physics_model import geometry
from app.services.benchmark_statistics import autocorrelation, split_rhat

STORE = ROOT / ".cache" / "layered-physics"
PAPER = ROOT / "manuscript_latest" / "review" / "overleaf"
EXECUTABLE = ROOT / "backend" / "build" / ("sse_layered.exe" if os.name == "nt" else "sse_layered")
PALETTE = {"ed": "#263B53", "q0": "#007C91", "qpi": "#D55E00", "ice": "#7A5195"}


def diagonal_observables(lx, ly, lz, basis_states=None):
    """Return diagonal observables for every sigma-z basis state."""
    area, n = lx * ly, lx * ly * lz
    plaquettes, pairs, signs = geometry(lx, ly)
    states = (2 * ((np.arange(1 << n)[:, None] >> np.arange(n)) & 1) - 1
              if basis_states is None else np.asarray(basis_states, dtype=int))
    pcount = len(plaquettes)
    total_plaquettes = pcount * lz
    ice = np.zeros(len(states)); rho4 = np.zeros(len(states)); plaq = np.zeros(len(states))
    pz = np.zeros((len(states), lz, 2))
    staggered = np.zeros((len(states), lz))
    for z in range(lz):
        layer = states[:, z*area:(z+1)*area]
        local = layer[:, plaquettes]
        px = local @ (np.array([1., -1., -1., 1.]) / 4.)
        py = local @ (np.array([1., 1., -1., -1.]) / 4.)
        sums = local.sum(axis=2)
        stable = (sums == 0) & (px*px + py*py > .5)
        ice += (sums == 0).sum(axis=1)
        rho4 += stable.sum(axis=1) - (pcount - stable.sum(axis=1))/3.
        plaq += np.prod(local, axis=2).sum(axis=1)
        pz[:, z, 0] = (px * signs).mean(axis=1)
        pz[:, z, 1] = (py * signs).mean(axis=1)
        pattern = np.fromfunction(lambda y, x: (-1.)**(x+y), (ly, lx)).reshape(-1)
        staggered[:, z] = layer @ pattern / area
    ice /= total_plaquettes; rho4 /= total_plaquettes; plaq /= total_plaquettes
    spin_q0 = staggered.mean(axis=1)
    spin_qpi = (staggered * ((-1.)**np.arange(lz))[None, :]).mean(axis=1)
    pmean0 = pz.mean(axis=1)
    pmeanpi = (pz * ((-1.)**np.arange(lz))[None, :, None]).mean(axis=1)
    dipole_q0 = (pmean0*pmean0).sum(axis=1)
    dipole_qpi = (pmeanpi*pmeanpi).sum(axis=1)
    return states, {
        "ice": ice, "rho4": rho4, "plaquette": plaq,
        "spin_q0_sq": spin_q0**2, "spin_qpi_sq": spin_qpi**2,
        "dipole_q0_sq": dipole_q0, "dipole_qpi_sq": dipole_qpi,
    }, (plaquettes, pairs, signs)


def exact_stack(params):
    """Complete finite-temperature ED for small periodic layered clusters."""
    lx, ly, lz = (int(params[k]) for k in ("Lx", "Ly", "Lz"))
    n = lx * ly * lz
    if n > 10:
        raise ValueError("The transparent dense ED reference is limited to 10 spins.")
    states, observables, (plaquettes, pairs, _) = diagonal_observables(lx, ly, lz)
    area = lx * ly
    diag = np.zeros(len(states), dtype=float)
    for z in range(lz):
        layer = states[:, z*area:(z+1)*area]
        local = layer[:, plaquettes]
        diag -= params["J0"] * np.prod(local, axis=2).sum(axis=1)
        diag += params["J1"] * (local[:, :, 0]*local[:, :, 2] + local[:, :, 1]*local[:, :, 3]).sum(axis=1)
        px = local @ (np.array([1., -1., -1., 1.]) / 4.)
        py = local @ (np.array([1., 1., -1., -1.]) / 4.)
        for a, b in pairs:
            diag += params["J2"] * (px[:, a]*px[:, b] + py[:, a]*py[:, b])
    for z in range(lz):
        diag += params["J3"] * np.sum(states[:, z*area:(z+1)*area] *
                    states[:, ((z+1) % lz)*area:((z+1) % lz+1)*area], axis=1)
    matrix = np.diag(diag)
    indices = np.arange(len(states))
    for site in range(n):
        matrix[indices, indices ^ (1 << site)] = -params["hx"]
    energies, vectors = eigh(matrix, driver="evd", check_finite=False)
    boltz = np.exp(-params["beta"] * (energies - energies[0]))
    weights = boltz / boltz.sum()
    probs = (vectors*vectors) @ weights
    result = {key: float(value @ probs) for key, value in observables.items()}
    result["E"] = float((energies @ weights) / n)
    result["cv"] = float(params["beta"]**2 * ((energies**2 @ weights) - (energies @ weights)**2) / n)
    residual = float(np.max(np.linalg.norm(matrix @ vectors-vectors*energies, axis=0)))
    return result, residual, energies


def _tau_and_error(streams):
    stats = [autocorrelation(x) for x in streams]
    taus = [v["tau"] for v in stats if v["tau"] is not None]
    tau = float(np.mean(taus)) if taus else None
    nmin = min(map(len, streams))
    block = max(1, min(max(32, math.ceil(10 * max(taus or [.5]))), max(1, nmin // 8)))
    means = np.asarray([float(x.mean()) for x in streams])
    weights = np.asarray([len(x) for x in streams], dtype=float); weights /= weights.sum()
    center = float(weights @ means)
    chain_se = math.sqrt(float(np.sum(weights**2 * (means-center)**2 / (1-weights)))) if len(streams)>1 else 0.
    blocks=[]; block_size=[]
    for x in streams:
        count=len(x)//block
        if count:
            blocks.extend(x[:count*block].reshape(count,block).mean(axis=1))
            block_size.extend([block]*count)
    if len(blocks)>1:
        # Equal-length delete-block jackknife. Chain-tail rows still enter the
        # point estimate, but not this variance component.
        all_rows=sum(map(len,streams)); full=sum(float(x.sum()) for x in streams)
        leave=[(full-value*size)/(all_rows-size) for value,size in zip(blocks,block_size)]
        jk=math.sqrt((len(leave)-1)/len(leave)*sum((x-np.mean(leave))**2 for x in leave))
    else: jk=0.
    acf_var=0.
    total=sum(map(len,streams))
    for x,s in zip(streams,stats):
        if s["tau"] is not None:
            acf_var += len(x)*float(np.var(x,ddof=1))*2*s["tau"]/total**2
    rhat=split_rhat(streams) if len(streams)>1 else None
    return center, max(chain_se,jk,math.sqrt(acf_var)), tau, int(total/block), rhat


def _run_chain(params, seed, algorithm, therm, measure, folder):
    folder.mkdir(parents=True, exist_ok=True)
    input_path, trace_path = folder / "input.in", folder / "trace.csv"
    fields={**params,"n_therm":therm,"n_measure":measure}
    if getattr(_run_chain,"reuse",False) and trace_path.exists():
        saved={line.split()[0]:float(line.split()[1]) for line in input_path.read_text(encoding="utf-8").splitlines() if len(line.split())==2}
        if any(not math.isclose(saved.get(k,math.nan),float(v),rel_tol=0,abs_tol=1e-12) for k,v in fields.items()):
            raise ValueError(f"Cannot reuse trace with different input parameters: {trace_path}")
    else:
        input_path.write_text("\n".join(f"{k} {v}" for k,v in fields.items())+"\n",encoding="utf-8")
        subprocess.run([str(EXECUTABLE),str(input_path),str(trace_path),algorithm,str(seed)],cwd=folder,check=True,capture_output=True,text=True)
    with trace_path.open(newline="",encoding="utf-8") as source:
        reader=csv.DictReader(source)
        rows=list(reader)
        fieldnames=reader.fieldnames
    return {key: np.asarray([float(row[key]) for row in rows]) for key in fieldnames}


def _style():
    plt.rcParams.update({
        "font.family":"serif", "font.serif":["STIX Two Text","STIXGeneral","DejaVu Serif"],
        "mathtext.fontset":"stix", "font.size":8, "axes.labelsize":8,
        "axes.titlesize":8, "legend.fontsize":7, "xtick.labelsize":7,
        "ytick.labelsize":7, "axes.linewidth":.7, "lines.linewidth":1.25,
        "lines.markersize":4, "pdf.fonttype":42, "ps.fonttype":42,
        "savefig.bbox":"tight", "savefig.pad_inches":.02,
    })


def run(args):
    if not EXECUTABLE.exists():
        raise RuntimeError("Build the trace engine first with `make -C backend/build layered`.")
    if args.lz < 2 or args.lz % 2 or args.lx % 2 or args.ly % 2 or args.lx*args.ly*args.lz>10:
        raise ValueError("Use even periodic Lx,Ly,Lz and at most 10 spins for this ED benchmark.")
    if args.chains < 4:
        raise ValueError("Use at least four independent chains for the reported uncertainty.")
    STORE.mkdir(parents=True,exist_ok=True)
    _run_chain.reuse=args.reuse
    run_id=f"{args.lx}x{args.ly}x{args.lz}-j0{args.j0:g}-b{args.beta:g}-j1{args.j1:g}-j2{args.j2:g}-h{args.hx:g}-{args.algorithm}"
    outdir=STORE/run_id
    outdir.mkdir(parents=True,exist_ok=True)
    j3_values=sorted(set(args.j3))
    records=[]; all_chains={}
    observables=("E","ice","rho4","plaquette","spin_q0_sq","spin_qpi_sq","dipole_q0_sq","dipole_qpi_sq")
    for point,j3 in enumerate(j3_values):
        params={"Lx":args.lx,"Ly":args.ly,"Lz":args.lz,"J0":args.j0,"J1":args.j1,"J2":args.j2,"J3":j3,"hx":args.hx,"beta":args.beta}
        ed, residual, energies=exact_stack(params)
        chain_data=[]
        for chain in range(args.chains):
            seed=args.seed + point*100003 + chain*7919
            traces=_run_chain(params,seed,args.algorithm,args.therm,args.measure,outdir/f"j3-{j3:+.8g}"/f"chain-{chain:02d}")
            chain_data.append(traces)
        all_chains[j3]=chain_data
        qmc={}
        for key in observables:
            mean,se,tau,blocks,rhat=_tau_and_error([c[key] for c in chain_data])
            qmc[key]={"mean":mean,"se":se,"tau_int_sweeps":tau,"blocks":blocks,"split_rhat":rhat}
        # The SSE factorial-moment energy variance is a separate estimator;
        # retain and report it, but do not use it as a tuning criterion.
        energy_series=[c["E"] for c in chain_data]
        pooled_e=float(np.mean([v for c in chain_data for v in c["E"]]))
        pooled_e2=float(np.mean([v for c in chain_data for v in c["E2"]]))
        scale=args.beta**2*args.lx*args.ly*args.lz
        cv_mean=scale*(pooled_e2-pooled_e**2)
        # Delta-method influence series for the nonlinear heat-capacity
        # estimator. The block/ACF/independent-chain uncertainty pipeline is
        # applied to this centered first-order influence function.
        cv_influence=[scale*(c["E2"]-2*pooled_e*c["E"]) for c in chain_data]
        _,cv_se,cv_tau,cv_blocks,cv_rhat=_tau_and_error(cv_influence)
        qmc["cv"]={"mean":cv_mean,"se":cv_se,"tau_int_sweeps":cv_tau,"blocks":cv_blocks,"split_rhat":cv_rhat}
        record={"J3":j3,"ED":ed,"QMC":qmc,"ed_residual":residual,"hilbert_dimension":len(energies)}
        records.append(record)
        print(f"J3={j3:+.4f}: E ED={ed['E']:.7f}, SSE={qmc['E']['mean']:.7f} +/- {qmc['E']['se']:.2g}; rho4={qmc['rho4']['mean']:.5f} +/- {qmc['rho4']['se']:.2g}",flush=True)

    # Exact operator mapping check on the complete finite spectrum.
    spectral={r["J3"]: exact_stack({"Lx":args.lx,"Ly":args.ly,"Lz":args.lz,"J0":args.j0,"J1":args.j1,"J2":args.j2,"J3":r["J3"],"hx":args.hx,"beta":args.beta})[2] for r in records}
    mapping_error=0.
    for value in set(abs(x) for x in j3_values):
        if value and value in spectral and -value in spectral:
            mapping_error=max(mapping_error,float(np.max(np.abs(spectral[value]-spectral[-value]))))
    pulls=[]
    for r in records:
        for key in (*observables,"cv"):
            ref=r["ED"].get(key); q=r["QMC"][key]
            if ref is not None and q["se"]:
                pulls.append({"J3":r["J3"],"observable":key,"z":(q["mean"]-ref)/q["se"]})
    source_paths=[ROOT/"scripts/layered_physics.py",ROOT/"backend/benchmarks/layered.cpp",ROOT/"backend/src/updates.cpp",ROOT/"backend/src/variables.cpp",ROOT/"backend/src/observables.cpp",ROOT/"backend/include/variables.hpp",ROOT/"backend/include/updates.hpp",ROOT/"backend/include/observables.hpp",ROOT/"backend/include/read_input.hpp",ROOT/"backend/app/services/physics_model.py",ROOT/"backend/app/services/benchmark_statistics.py"]
    source_hashes={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}
    revision=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip()
    report={"schema":"layered-sign-benchmark-v1","geometry":[args.lx,args.ly,args.lz],"boundary":"periodic in x,y,z; even Lz",
            "parameters":{"J0":args.j0,"J1":args.j1,"J2":args.j2,"hx":args.hx,"beta":args.beta,"J3_values":j3_values},
            "sampling":{"algorithm":args.algorithm,"chains":args.chains,"thermalization_sweeps":args.therm,"measurement_sweeps_per_chain":args.measure,"seed_base":args.seed,
                        "error":"maximum of block jackknife, autocorrelation variance, independent-chain HC2 SE","tau_convention":"tau_int=1/2+sum rho(k), units of complete sweeps"},
            "ed_sign_mapping_max_spectral_difference":mapping_error,
            "maximum_abs_qmc_ed_residual_in_reported_se":max((abs(v["z"]) for v in pulls),default=None),
            "all_observable_pulls_within_4_se":all(abs(v["z"])<=4 for v in pulls),"pulls":pulls,
            "provenance":{"git_revision":revision,"working_tree_dirty":bool(subprocess.check_output(["git","status","--porcelain"],cwd=ROOT,text=True).strip()),
                          "binary_sha256":hashlib.sha256(EXECUTABLE.read_bytes()).hexdigest(),"source_sha256":source_hashes,
                          "python":platform.python_version(),"numpy":np.__version__},"runs":records}
    (outdir/"summary.json").write_text(json.dumps(report,indent=2,allow_nan=False),encoding="utf-8")
    with (outdir/"observables.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.writer(f); writer.writerow(["J3","observable","ED","SSE","SE","tau_int_sweeps","split_rhat","z_score"])
        for row in records:
            for key in (*observables,"cv"):
                q=row["QMC"][key]; ed=row["ED"].get(key)
                z=(q["mean"]-ed)/q["se"] if ed is not None and q["se"] else None
                writer.writerow([row["J3"],key,ed,q["mean"],q["se"],q["tau_int_sweeps"],q["split_rhat"],z])
    make_figure(records,outdir,args)
    make_all_observables_figure(records,outdir,args)
    write_table(records,PAPER/"tables"/"layered-sign-benchmark.tex",args)
    (PAPER/"figures").mkdir(parents=True,exist_ok=True)
    for stem in ("layered-physics-results","layered-all-observables"):
        for suffix in ("pdf","png","svg"):
            source=outdir/f"{stem}.{suffix}"
            if source.exists(): (PAPER/"figures"/source.name).write_bytes(source.read_bytes())
    print(f"ED spectral sign-mapping max |delta E_n| = {mapping_error:.3e}")
    print(f"Maximum |ED-SSE| / reported SE = {report['maximum_abs_qmc_ed_residual_in_reported_se']:.2f}")
    print(f"Local results: {outdir}")


def make_figure(records,outdir,args):
    _style()
    fig,axs=plt.subplots(2,2,figsize=(7.15,4.75),sharex=True)
    x=np.asarray([r["J3"] for r in records])
    ed_color="#263B53"; qmc_color="#007C91"
    common=[("E",r"(a) Energy per pseudospin",r"$\langle H\rangle/N_s$"),
            ("rho4",r"(b) Polarized-ice diagnostic",r"$\langle\rho_4\rangle$")]
    for ax,(key,title,ylabel) in zip(axs[0],common):
        edvals=[r["ED"][key] for r in records]
        qvals=[r["QMC"][key]["mean"] for r in records]
        errs=[r["QMC"][key]["se"] for r in records]
        ax.plot(x,edvals,color=ed_color,marker="s",ms=4.,lw=1.1,label="ED, exact trace",zorder=2)
        ax.errorbar(x,qvals,yerr=errs,fmt="o",ms=5.,mfc="white",mec=qmc_color,
                    ecolor=qmc_color,capsize=2.2,lw=.85,label="SSE QMC, 1 SE",zorder=4)
        ax.set_title(title,loc="left",pad=5)
        ax.set_ylabel(ylabel)

    channel_specs=[
        (axs[1,0],r"(c) Dipole stacking channels",("dipole_q0_sq","dipole_qpi_sq"),r"$S_P(q_z)$"),
        (axs[1,1],r"(d) Staggered-spin stacking channels",("spin_q0_sq","spin_qpi_sq"),r"$S_m(q_z)$"),
    ]
    for ax,title,keys,ylabel in channel_specs:
        for key,color,marker,label in ((keys[0],PALETTE["q0"],"o",r"$q_z=0$"),
                                       (keys[1],PALETTE["qpi"],"^",r"$q_z=\pi$")):
            edvals=[r["ED"][key] for r in records]
            qvals=[r["QMC"][key]["mean"] for r in records]
            errs=[r["QMC"][key]["se"] for r in records]
            ax.plot(x,edvals,color=color,marker="s",ms=3.7,lw=1.0,label=f"ED {label}",zorder=2)
            ax.errorbar(x,qvals,yerr=errs,fmt=marker,ms=4.2,mfc="white",mec=color,
                        ecolor=color,capsize=2.0,lw=.8,label=f"SSE {label}",zorder=4)
        ax.set_title(title,loc="left",pad=5)
        ax.set_ylabel(ylabel)
    for ax in axs.flat:
        ax.set_xlabel(r"Interlayer coupling $J_3/J_0$")
        ax.set_xticks(x)
        ax.tick_params(direction="in",top=True,right=True,length=3.0,width=.65,labelsize=7.2)
        ax.margins(x=.08,y=.12)
        for spine in ax.spines.values(): spine.set_linewidth(.65)
    from matplotlib.lines import Line2D
    handles=[
        Line2D([0],[0],color=ed_color,marker="s",lw=1.0,ms=4,label="ED, exact trace"),
        Line2D([0],[0],color=qmc_color,marker="o",mfc="white",mec=qmc_color,
               lw=0,ms=4.5,label="SSE QMC (bars: 1 SE)"),
        Line2D([0],[0],color=PALETTE["q0"],lw=1.3,label=r"$q_z=0$ channel"),
        Line2D([0],[0],color=PALETTE["qpi"],lw=1.3,label=r"$q_z=\pi$ channel"),
    ]
    fig.legend(handles=handles,frameon=False,ncol=4,loc="upper center",
               bbox_to_anchor=(.5,.955),handlelength=1.7,columnspacing=1.4,
               fontsize=6.9)
    title=(r"Periodic $%dx%dx%d$ stack; $\beta J_0=%g$, $J_1/J_0=%g$, $J_2/J_0=%g$, $h_x/J_0=%g$"
           % (args.lx,args.ly,args.lz,args.beta,args.j1,args.j2,args.hx))
    fig.suptitle(title,fontsize=8.1,y=1.005)
    fig.subplots_adjust(left=.10,right=.99,bottom=.11,top=.84,hspace=.32,wspace=.25)
    for suffix in ("pdf","svg","png"):
        fig.savefig(outdir/f"layered-physics-results.{suffix}",dpi=600 if suffix=="png" else None)
    plt.close(fig)


def write_table(records,path,args):
    path.parent.mkdir(parents=True,exist_ok=True)
    rows=[]
    for r in records:
        e=r["ED"];q=r["QMC"]
        rows.append("{:+.2f} & {:.6f} & {:.6f} $\\pm$ {:.6f} & {:.6f} & {:.6f} $\\pm$ {:.6f} & {:.6f} & {:.6f} $\\pm$ {:.6f} & {:.6f} & {:.6f} $\\pm$ {:.6f} \\\\".format(
            r["J3"],e["E"],q["E"]["mean"],q["E"]["se"],e["rho4"],q["rho4"]["mean"],q["rho4"]["se"],
            e["dipole_q0_sq"],q["dipole_q0_sq"]["mean"],q["dipole_q0_sq"]["se"],e["dipole_qpi_sq"],q["dipole_qpi_sq"]["mean"],q["dipole_qpi_sq"]["se"]))
    caption=(f"Complete-ED and SSE estimates for the periodic ${args.lx}\\times{args.ly}\\times{args.lz}$ stack at "
             f"$\\beta J_0={args.beta:g}$, $J_1/J_0={args.j1:g}$, $J_2/J_0={args.j2:g}$, and $h_x/J_0={args.hx:g}$. "
             "ED has no sampling uncertainty; SSE values are followed by one correlation-aware standard error. "
             "$S_P(0)$ and $S_P(\\pi)$ are the squared layer-averaged dipole amplitudes at the indicated stacking wave vectors.")
    content=("\\begin{table*}[t]\n\\centering\\scriptsize\n"
        f"\\caption{{{caption}}}\n"
        "\\label{tab:layered-sign-benchmark}\n\\resizebox{\\textwidth}{!}{%\n\\begin{tabular}{r r r r r r r r r}\n\\hline\\hline\n"
        "$J_3/J_0$ & $e_{\\rm ED}$ & $e_{\\rm SSE}$ & $\\rho_{4,\\rm ED}$ & $\\rho_{4,\\rm SSE}$ & $S_P(0)_{\\rm ED}$ & $S_P(0)_{\\rm SSE}$ & $S_P(\\pi)_{\\rm ED}$ & $S_P(\\pi)_{\\rm SSE}$ \\\\\n\\hline\n"+"\n".join(rows)+"\n\\hline\\hline\n\\end{tabular}}\n\\end{table*}\n")
    path.write_text(content,encoding="utf-8")


def make_all_observables_figure(records,outdir,args):
    _style()
    panels=[("E",r"Energy per pseudospin",r"$\langle H\rangle/N_s$"),
            ("cv",r"Specific heat (low precision)",r"$C/N_s$"),
            ("ice",r"Two-up/two-down fraction",r"$f_{2:2}$"),
            ("rho4",r"Polarized-ice diagnostic",r"$\langle\rho_4\rangle$"),
            ("plaquette",r"Plaquette product",r"$\langle A_{\rm plaq}\rangle$"),
            ("spin_q0_sq",r"Staggered spin, $q_z=0$",r"$S_m(0)$"),
            ("spin_qpi_sq",r"Staggered spin, $q_z=\pi$",r"$S_m(\pi)$"),
            ("dipole_q0_sq",r"Dipole, $q_z=0$",r"$S_P(0)$"),
            ("dipole_qpi_sq",r"Dipole, $q_z=\pi$",r"$S_P(\pi)$")]
    fig,axs=plt.subplots(3,3,figsize=(7.15,6.05),sharex=True,constrained_layout=True)
    x=np.asarray([r["J3"] for r in records])
    for i,(ax,(key,title,ylabel)) in enumerate(zip(axs.flat,panels)):
        ed=np.asarray([r["ED"][key] if key!="cv" else r["ED"]["cv"] for r in records])
        q=np.asarray([r["QMC"][key]["mean"] for r in records])
        se=np.asarray([r["QMC"][key]["se"] for r in records])
        ax.plot(x,ed,color=PALETTE["ed"],marker="s",ms=3.1,lw=1.,label="ED, exact trace")
        ax.errorbar(x,q,yerr=se,fmt="o",ms=3.,mfc="white",mec=PALETTE["q0"],
                    ecolor=PALETTE["q0"],capsize=1.6,lw=.75,label="SSE, 1 SE")
        ax.set_title(f"({chr(97+i)}) {title}",loc="left",pad=3);ax.set_ylabel(ylabel)
        ax.tick_params(direction="in",top=True,right=True,length=2.7,width=.6)
        ax.margins(x=.08)
        if i//3==2: ax.set_xlabel(r"$J_3/J_0$")
    handles,labels=axs.flat[0].get_legend_handles_labels()
    fig.legend(handles,labels,frameon=False,loc="upper center",ncol=2,bbox_to_anchor=(.5,1.025))
    fig.suptitle(fr"$ {args.lx}\times{args.ly}\times{args.lz}$; $\beta J_0={args.beta:g}$; $J_1/J_0={args.j1:g}$; $J_2/J_0={args.j2:g}$; $h_x/J_0={args.hx:g}$",fontsize=8,y=1.065)
    for suffix in ("pdf","svg","png"):
        fig.savefig(outdir/f"layered-all-observables.{suffix}",dpi=600 if suffix=="png" else None)
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--lx",type=int,default=2);p.add_argument("--ly",type=int,default=2);p.add_argument("--lz",type=int,default=2)
    p.add_argument("--j0",type=float,default=1.);p.add_argument("--j1",type=float,default=.3);p.add_argument("--j2",type=float,default=.025);p.add_argument("--hx",type=float,default=.2)
    p.add_argument("--beta",type=float,default=4.);p.add_argument("--j3",type=float,nargs="+",default=[-.3,-.2,0.,.2,.3])
    p.add_argument("--chains",type=int,default=8);p.add_argument("--therm",type=int,default=5000);p.add_argument("--measure",type=int,default=25000)
    p.add_argument("--algorithm",choices=("vertex","premarked","canonical","symmetry"),default="symmetry");p.add_argument("--seed",type=int,default=20261018)
    p.add_argument("--reuse",action="store_true",help="Analyze matching existing local traces without rerunning SSE.")
    args=p.parse_args();run(args)


if __name__=="__main__": main()
