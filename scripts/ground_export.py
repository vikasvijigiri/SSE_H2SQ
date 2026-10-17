# Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
"""Export a completed frozen ground study to an ignored Overleaf folder."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import statistics
import sys
import urllib.request

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'backend'))


def number(value):
    return '--' if value is None else f'{value:.6g}'


def tex_number(value):
    text=number(value)
    if 'e' not in text:return text
    mantissa,exponent=text.split('e')
    return mantissa+r'\times10^{'+str(int(exponent))+'}'


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('identity')
    parser.add_argument('--name',default='ground-state-validation');parser.add_argument('--reuse-figures',action='store_true');args=parser.parse_args()
    if not re.fullmatch(r'[A-Za-z0-9_-]+',args.name):raise SystemExit('Use a simple folder name.')
    with urllib.request.urlopen(f'http://127.0.0.1:8000/api/ground-validation/{args.identity}',timeout=60) as response:job=json.load(response)
    if job['status']!='complete':raise SystemExit('Complete the frozen study before exporting final evidence.')
    from app.services.ground_validation import csv_result,figure_result,PHYSICS_KEYS,STORE
    folder=ROOT/'docs/manuscript/benchmark'/args.name
    existing=folder/'metadata/experiment.json'
    if existing.exists() and json.loads(existing.read_text())['id']!=job['id']:
        raise SystemExit('This folder holds a different study. Choose another --name; original results are preserved.')
    for name in ['figures','tables','metadata','sections','bibliography']:(folder/name).mkdir(parents=True,exist_ok=True)
    (folder/'metadata/experiment.json').write_text(json.dumps(job,indent=2,allow_nan=False),encoding='utf-8')
    (folder/'tables/observables.csv').write_text(csv_result(job),encoding='utf-8')
    rows=[];panels=0;figure_sections=[]
    for point in job['points']:
        p=point['params'];levels=point['levels'];last=levels[-1] if levels else None
        histogram=last['methods'].get('projector',{}).get('histogram',{}) if last else {}
        rows.append([point['id'],point['kind'],p['Lx'],p['Ly'],p['J0'],p['J1'],p['J2'],p['hx'],last['projection'] if last else None,point['verdict'],histogram.get('fidelity_ground'),histogram.get('tvd_ground'),histogram.get('missing_ground_weight'),'; '.join(point.get('reasons',[])) or point.get('reference_error','')])
        if not levels or not any(l['methods'] for l in levels):continue
        choices=[(key,'mean') for key in PHYSICS_KEYS+['E_count','E2','E4']]
        choices += [(key,quantity) for key in ['E','ct0','stripe2'] for quantity in ['sampling_difference','tau_int','precision']]
        if not point['kind'].startswith('classical'):choices += [('E','fidelity')]
        if histogram:choices += [('E',key) for key in ['state_infidelity','state_tvd','state_missing']]
        for key,quantity in choices:
            for form in ['pdf','svg','png']:
                target=folder/f'figures/{point["id"]}-{key}-{quantity}.{form}'
                if not args.reuse_figures or not target.is_file():target.write_bytes(figure_result(job,point['id'],key,quantity,form))
            panels+=1
        title=point['label'].replace('_',r'\_')
        figure_sections.append(r'\subsection{'+title+'}\n'+r'\noindent '+f'{p["Lx"]} by {p["Ly"]}; '+r'$J_1='+number(p['J1'])+r', J_2='+number(p['J2'])+r', h='+number(p['hx'])+r'$. Verdict: \textbf{'+point['verdict']+'}.\n')
        assets=[point['id']+'-'+key+'-mean.pdf' for key in ['E','mx','stripe2','ct0']]
        if histogram:
            assets += [point['id']+'-E-'+quantity+'.pdf' for quantity in ['state_infidelity','state_tvd']]
        figure_sections.append(r'\par\smallskip'+'\n')
        for i in range(0,len(assets),2):
            figure_sections.append(r'\noindent '+r'\hfill'.join(r'\includegraphics[width=.48\linewidth]{figures/'+name+'}' for name in assets[i:i+2])+r'\par\smallskip'+'\n')
        figure_sections.append(r'\clearpage'+'\n')
        print(f'Exported {point["id"]}: {panels} panels',flush=True)
    with (folder/'tables/cases.csv').open('w',newline='',encoding='utf-8') as output:
        writer=csv.writer(output);writer.writerow(['point','kind','Lx','Ly','J0','J1','J2','hx','T_or_beta','verdict','QMC_ground_fidelity','TVD_ground','missing_ground_weight','reasons']);writer.writerows(rows)
    with (folder/'tables/averaging.csv').open('w',newline='',encoding='utf-8') as output:
        writer=csv.writer(output);writer.writerow(['point','method','runs','rows','observable','SE','measured_ratio','expected_ratio'])
        for point in job['points']:
            if not point['levels']:continue
            for method,summary in point['levels'][-1]['methods'].items():
                for history in summary.get('precision_history',[]):
                    for key,v in history['observables'].items():writer.writerow([point['id'],method,history['runs'],history['rows'],key,v['error'],v['SE_ratio'],history['expected_SE_ratio']])
    table=[r'\begin{longtable}{lrrrrrrl}',r'Case & $N$ & $J_1$ & $J_2$ & $h$ & $T/\beta$ & $F_{\rm QMC,0}$ & Verdict \\ \hline']
    for row in rows:table.append(' & '.join([row[0],str(row[2]*row[3]),*[number(row[i]) for i in [5,6,7,8,10]],row[9]])+r' \\')
    table.append(r'\end{longtable}');(folder/'tables/cases.tex').write_text('\n'.join(table),encoding='utf-8')
    estimates=[r'{\footnotesize\begin{longtable}{llrrl}',r'Case / method & Quantity & QMC $\pm$ SE & Ground / exact & Status \\ \hline\endhead']
    for point in job['points']:
        if not point['levels']:continue
        for method,summary in point['levels'][-1]['methods'].items():
            for key in summary['observables']:
                value=summary['observables'][key]
                estimate='$'+tex_number(value['mean'])+(r'\pm'+tex_number(value['error']) if value['error'] is not None else r'\;\mathrm{(SE\ unresolved)}')+'$'
                reference='$'+tex_number(value['ground_ed'])+'$' if value['ground_ed'] is not None else '--'
                estimates.append(' & '.join([point['id']+'/'+method,key.replace('_',r'\_'),estimate,reference,value['status']])+r' \\')
    estimates.append(r'\end{longtable}}');(folder/'tables/final-estimates.tex').write_text('\n'.join(estimates),encoding='utf-8')
    (folder/'sections/figures.tex').write_text(''.join(figure_sections),encoding='utf-8')
    (folder/'sections/methods.tex').write_text(r'''\section{Hamiltonian and ground-state target}
The periodic even checkerboard lattice uses Pauli variables $\sigma_i^z=\pm1$. The cyclic black-plaquette order is
\[ (0,1,2,3)=(LL,UL,UR,LR). \]
\begin{equation}
 H=-J_0\sum_p\sigma_{p0}\sigma_{p1}\sigma_{p2}\sigma_{p3}
 +J_1\sum_p(\sigma_{p0}\sigma_{p2}+\sigma_{p1}\sigma_{p3})
 +J_2\sum_{\langle AB\rangle}\widetilde{\mathbf P}_A\cdot\widetilde{\mathbf P}_B
 -h\sum_i\sigma_i^x .
\end{equation}
The tilde dipoles omit the molecular B sign. Thus the displayed $+J_2$ term equals the manuscript's $-J_2\mathbf P_A\cdot\mathbf P_B$ with $\mathbf P_B=-\widetilde{\mathbf P}_B$. Directional periodic multiplicities are retained even when endpoint pairs coincide in narrow cells. Spin moments use $S^z=\sigma^z/2$; transverse magnetization uses Pauli $\sigma^x$ per spin. Ground response functions are outside this report's scope.

We sample $\langle+|e^{-TH}|+\rangle$ with free imaginary-time endpoints, corresponding to the normalized pure state $e^{-TH/2}|+\rangle$. This is a projector, not a thermal trace. Uniform nonzero $h$ is gauged to $|h|$; the physical transverse estimator retains its sign. This product trial and open boundaries follow the established projector construction~\cite{sandvik2019}.

The 12 required physical static quantities are energy, staggered second/fourth moments, transverse magnetization, plaquette product, ice-rule fraction, $\rho$, stripe second/fourth moments, polarization second/fourth moments, and the equal-time dipole quantity $c_{t0}$. Their definitions and normalizations are retained with the frozen study. Count energy, second/fourth energy moments and four raw estimators provide additional audits.

For a sampled $z$ configuration, define
\begin{align}
 m_z&=\frac{1}{2N}\sum_{x,y}(-1)^{x+y}\sigma_{x,y}, &
 s_2&=\left(\frac{\sum_{x,y}(-1)^x\sigma_{x,y}}{2N}\right)^2
       +\left(\frac{\sum_{x,y}(-1)^y\sigma_{x,y}}{2N}\right)^2,\\
 \widetilde P_{p,x}&=\tfrac14(\sigma_0-\sigma_1-\sigma_2+\sigma_3), &
 \widetilde P_{p,y}&=\tfrac14(\sigma_0+\sigma_1-\sigma_2-\sigma_3).
\end{align}
Staggered moments are $\langle m_z^2\rangle,\langle m_z^4\rangle$; stripe moments are $\langle s_2\rangle,\langle s_2^2\rangle$. Physical dipoles are $\mathbf P_p=\eta_p\widetilde{\mathbf P}_p$ with $\eta_p=+1$ on A and $-1$ on B. Polarization moments use $p_2=|N_p^{-1}\sum_p\mathbf P_p|^2$ and $p_2^2$. The plaquette product is the mean four-spin product; $f_{\rm ice}$ is the fraction with $\sum_{i\in p}\sigma_i=0$ and $\rho=(4f_{\rm ice}-1)/3$. The measured $c_{t0}=N_p^{-1}\sum_p\langle\widetilde P_{p,x}^2\rangle$ is local, equal-time and x-component only. The transverse mean is $N^{-1}\sum_i\langle\sigma_i^x\rangle$. Energy moments are $\langle(H/N)^k\rangle$ for $k=1,2,4$.

\section{Estimators and energy shifts}
The central insertion index has exact conditional law $p\sim\mathrm{Binomial}(n,1/2)$. Diagonal estimators integrate constant propagation segments with these probabilities. The central off-diagonal event intensity is $(n/T)\mathrm{Binomial}(n-1,r-1;1/2)$ at event rank $r$, not the boundary-averaged total off-event count. At a fixed cluster graph, distinct-component field endpoints give conditional off probability $1/2$; same-component endpoints retain their actual off status, including frustrated contributions. Component-orientation averages give $\langle m^2\rangle_G=\sum_c a_c^2$ and $\langle m^4\rangle_G=3(\sum_c a_c^2)^2-2\sum_c a_c^4$. Model-specific formulas were checked against exhaustive component enumeration, following the improved-estimator principle~\cite{sandvik2011}. Raw estimator audits are retained.

Write $K=C-H$. The physical count energy is $(C-\langle n\rangle/T)/N$. Here $C=N_p|J_0|+2N_p|J_1|+\sum_b|J_2J_b|+J_2c_{\rm geom}+N|h|$. The geometric constant is physical and is included exactly once. Higher physical energy moments use falling factorials and the binomial expansion of $(C-K)^k/N^k$. Midpoint $\langle H_z\rangle/N-hm_x$ independently audits the count energy. A fixed model-norm cutoff and Poisson tail provide explicit absolute observable-bias bounds; zero observed saturations alone is not a proof of negligible truncation.

\section{Balance, support and convergence}
Insertion/removal uses the full padded-string Metropolis ratio; local identity probability is positive. Spatial premarking and canonical updates are motivated by frustrated quantum-cluster constructions~\cite{biswas2019}; all rules are checked against this Hamiltonian's own vertex weights. Cluster component orientations are exact conditional heat baths: pair and plaquette products, graph motif weights and uniform boundary amplitudes are preserved. Whole-worldline canonical proposals use all propagated diagonal factors. Independent even/odd spin flips preserve the actual Hamiltonian and all vertex weights. The ordered composition is stationary; it need not itself be reversible.

For fixed finite cutoff, a positive-probability path removes all interaction diagonal vertices while other steps idle. With only field vertices, open time segments can be oriented to make all field events diagonal, which can then be removed. Empty strings permit arbitrary boundary-spin changes. Reversing these elementary positive-probability steps connects every positive-weight configuration. Identity probabilities give aperiodicity. This establishes finite-state support under the stated model and boundaries, not a useful uniform mixing-time bound. Near small gaps or first-order coexistence, convergence can still be arbitrarily slow.

\section{Independent references and statistical criteria}
ED includes the complete symmetry-invariant Hilbert subspace containing the unique positive Perron--Frobenius ground state for nonzero field. Every group action is checked on the actual Hamiltonian; full-space ground and all reduced eigenpair residuals are checked. Finite-projection references use the complete spectrum of that subspace, never an assumed thermal low-state truncation. All sampled central configuration probabilities are averaged over known exact geometry orbits. Their positive-amplitude reconstruction is compared with the exact projected and ground states using Bhattacharyya fidelity, total variation and missing reference weight, with no ED-informed pseudocounts. Whole-independent-run bootstrap errors are conditional on observed support and are descriptive. Reconstructed variational energies give another diagnostic, not a finite-sample fidelity certificate.

Independent classical exchange Monte Carlo~\cite{hukushima1996} uses separately implemented geometry and Metropolis energies and is compared with thermal SSE and complete classical enumeration at $h=0$. Replica round trips are recorded. Such checks do not validate a nonzero-field quantum ground state. Larger product-state limits similarly validate estimators and size handling only.

Reported bars are one standard error. The SE is the maximum of blocked, per-run autocorrelation, and independent-run HC2 estimates. Rank/folded split $\widehat R$, ESS, block adequacy, complete-sweep autocorrelation and per-stage wall timings are retained. ESS/runtime divides by summed measured run wall times, including thermalization, and is not a CPU-time measurement. Autocorrelation is measured for the complete update composition, not isolated update kernels. Prefix four/eight-run comparisons retain SE increases. The frozen rule requires finite-ensemble agreement within four SE, requested precision, controlled projection bias, two long projection lengths, and state-distribution thresholds. These are review criteria, not simultaneous 100\% confidence bounds. Projection fidelity predicted by ED and sampled QMC-state fidelity are separate plots.
''',encoding='utf-8')
    counts={verdict:sum(p['verdict']==verdict for p in job['points']) for verdict in ['GO','NO GO']}
    summaries=[s for p in job['points'] for level in p['levels'] for s in level['methods'].values()]
    observed=[v for s in summaries for v in s['observables'].values()]
    stochastic=[v for v in observed if v['tau_int'] is not None]
    ratios=[v['SE_ratio'] for p in job['points'] if p['levels'] for s in p['levels'][-1]['methods'].values() for history in s.get('precision_history',[])[-1:] for v in history['observables'].values() if v['SE_ratio'] is not None and v['tau_int'] is not None]
    random=[p for p in job['points'] if p['label'].startswith('Random')]
    evidence=dict(rows=sum(s['retained_samples'] for s in summaries),cases=counts,
        random_passed=sum(p['verdict']=='GO' for p in random),random_total=len(random),
        maximum_stochastic_ED_pull=max((abs(v['pull']) for v in stochastic if v['pull'] is not None),default=None),
        maximum_rhat=max((v['rhat'] for v in stochastic if v['rhat'] is not None),default=None),
        minimum_ess=min((v['ess'] for v in stochastic),default=None),
        maximum_tau=max((v['tau_int'] for v in stochastic),default=None),
        stochastic_SE_comparisons=len(ratios),stochastic_SE_decreased=sum(r<1 for r in ratios),
        stochastic_SE_ratio_median=statistics.median(ratios) if ratios else None)
    (folder/'metadata/validation-summary.json').write_text(json.dumps(evidence,indent=2,allow_nan=False),encoding='utf-8')
    results=r'\section{Frozen study and verdicts}'+'\n'+f'Study {job["id"]}; revision '+r'\texttt{'+job['provenance']['revision'][:12]+'}. '+f'{job["completed_runs"]} independent runs and {evidence["rows"]:,} measured rows completed. '+r'Universal verdict: \textbf{NO GO}. Requested random-domain verdict: \textbf{'+job['verdict']['random_ground_domain']+'}. '+f'{counts["GO"]} cases pass their stated scope; {counts["NO GO"]} do not. '+r'Failed, stress and unsupported cases remain in the tables and figures. All individual quantities, diagnostics and reasons are in the accompanying CSV and JSON. These data support only the reported finite-cell claims.'+'\n\n'
    if random:results+=f'{evidence["random_passed"]} of {evidence["random_total"]} simultaneous random-coupling cases pass all required checks. '
    results+=f'The maximum finite-reference standardized residual among nonconstant estimators is {number(evidence["maximum_stochastic_ED_pull"])} SE. '+r'This alone is insufficient: sampling diagnostics fail even when mean residuals fit broad error bars. Deterministic conditional estimators have only floating-point residuals; dividing those by roundoff-scale SE is not a meaningful stochastic test. '
    results+=f'The largest measured autocorrelation is {number(evidence["maximum_tau"])} full sweeps, minimum ESS {number(evidence["minimum_ess"])}, and maximum '+r'$\widehat R$ '+number(evidence['maximum_rhat'])+'. Long autocorrelation estimates in flagged runs remain unresolved.\n\n'
    if ratios:results+=f'At the final projection/equilibrium length, SE decreases in {sum(r<1 for r in ratios)} of {len(ratios)} nonconstant four/eight-run comparisons; the median ratio is {statistics.median(ratios):.4f}. '+r'The ideal independent-run ratio is $1/\sqrt{2}\simeq0.7071$. Increases are retained. These comparisons do not demonstrate convergence of a flagged slow mode.'+'\n\n'
    results+=r'The case precision gate uses the 12 physical static quantities listed in the methods. Count energy, higher energy moments and four raw estimators are additional audits; inefficient raw-estimator precision does not veto a converged conditional estimate. Analytic product-state controls and classical checks have their own stated scope.'+'\n'+r'\input{tables/cases.tex}'+'\n'
    distributions=[(p,p['levels'][-1]) for p in job['points'] if p['levels'] and p['levels'][-1]['methods'].get('projector',{}).get('histogram')]
    if distributions:
        p,level=min(distributions,key=lambda pair:pair[1]['methods']['projector']['histogram']['fidelity_ground'])
        summary=level['methods']['projector'];h=summary['histogram'];e=summary['observables']['E']
        results+=r'\subsection{Why energy agreement is insufficient}'+'\n'+f'The lowest sampled ground-state fidelity occurs in case {p["id"]} at total projection {level["projection"]:g}. '+r'Energy per spin is $'+tex_number(e['mean'])+r'\pm'+tex_number(e['error'])+r'$; finite-projector ED gives $'+tex_number(e['projected_ed'])+r'$ and ground ED gives $'+tex_number(e['ground_ed'])+r'$. Sampled-state fidelity is $'+tex_number(h['fidelity_projected'])+r'$ to the finite projector, but $'+tex_number(h['fidelity_ground'])+r'$ to the ground state. Ground-distribution total variation is $'+tex_number(h['tvd_ground'])+r'$. Correct sampling of a finite projector does not imply sufficient ground projection.'+'\n\n'
        slow=[(p,l,s) for p,l in distributions for s in [l['methods']['projector']] if s['observables']['mz2']['tau_int'] is not None]
        if slow:
            p,l,s=max(slow,key=lambda item:item[2]['observables']['mz2']['tau_int']);v=s['observables']['mz2'];e=s['observables']['E']
            results+=f'The slowest final-level staggered mode is case {p["id"]}: '+r'$\tau_{\rm int}='+number(v['tau_int'])+r'$, ESS $'+number(v['ess'])+r'$, and $\widehat R='+number(v['rhat'])+r'$. Its energy residual is only $'+number(e['pull'])+r'$ SE from finite-projector ED. This mode requires improved mixing; broad error bars cannot justify GO.'+'\n'
    (folder/'sections/results.tex').write_text(results,encoding='utf-8')
    (folder/'sections/limitations.tex').write_text(r'''\section{Limits and next algorithms}
No finite benchmark proves correctness for every real coupling and every lattice size with 100\% confidence. Interacting ED is limited here to at most 16 spins; even periodic 2D geometry, uniform field, $J_3=0$ and the stated precision define this study. No thermodynamic phase boundary, real-time dynamics, or quantum ground-state susceptibility is certified. $h=0$ may have a degenerate ground manifold and needs explicit manifold/ensemble treatment. Very small nonzero fields or gaps can exceed any fixed projection or sampling budget; those cases must remain NO GO.

Ground projection with a symmetry-compatible trial can suppress irrelevant charge-sector contamination of a low-temperature thermal trace. It does not remove spatial free-energy barriers. Candidate improvements include projection/coupling replica exchange, extended-cell/directed-loop updates tailored to this plaquette model, and exact Metropolis-corrected self-learning proposals~\cite{sandvik2019,liu2017}. These are proposals, not results of this study. Learned proposals must preserve support, include the correct proposal/target ratio and freeze adaptation before production. Evaluate gains through independent replicas, tunnelling, ESS per CPU second and distribution agreement, not energy alone. An attempted bridge-only transverse estimator failed the frustrated small-cell exact test and was excluded; current estimators retain the missing within-component frustrated events.
''',encoding='utf-8')
    (folder/'bibliography/references.bib').write_text(r'''@article{sandvik2011,author={Sandvik, Anders W.},title={Computational Studies of Quantum Spin Systems},journal={AIP Conference Proceedings},volume={1297},pages={135--338},year={2010},doi={10.1063/1.3518900},eprint={1101.3281},archivePrefix={arXiv}}
@incollection{sandvik2019,author={Sandvik, Anders W.},title={Stochastic Series Expansion Methods},booktitle={Many-Body Methods for Real Materials},year={2019},eprint={1909.10591},archivePrefix={arXiv}}
@misc{biswas2019,author={Biswas, Sounak and Damle, Kedar},title={Efficient quantum cluster algorithms for frustrated transverse field Ising antiferromagnets and Ising gauge theories},year={2019},eprint={1812.05326},archivePrefix={arXiv}}
@article{hukushima1996,author={Hukushima, Koji and Nemoto, Koji},title={Exchange Monte Carlo Method and Application to Spin Glass Simulations},journal={Journal of the Physical Society of Japan},volume={65},pages={1604--1608},year={1996},doi={10.1143/JPSJ.65.1604},eprint={cond-mat/9512035},archivePrefix={arXiv}}
@article{liu2017,author={Liu, Junwei and Qi, Yang and Meng, Zi Yang and Fu, Liang},title={Self-Learning Monte Carlo Method},journal={Physical Review B},volume={95},pages={041101},year={2017},doi={10.1103/PhysRevB.95.041101},eprint={1610.03137},archivePrefix={arXiv}}
''',encoding='utf-8')
    (folder/'main.tex').write_text(r'''\documentclass[11pt]{article}
\usepackage[margin=22mm]{geometry}
\usepackage{amsmath,amssymb,graphicx,longtable,booktabs,hyperref}
\setlength{\emergencystretch}{2em}
\title{Finite-lattice ground-state validation of checkerboard projector QMC}
\author{Vikas Vijigiri}
\date{\today}
\begin{document}\maketitle
\noindent Technical validation report. Scope and failed cases are explicit; this report does not claim universal or 100\% validation.
\input{sections/methods.tex}
\input{sections/results.tex}
\input{sections/limitations.tex}
\bibliographystyle{unsrt}\bibliography{bibliography/references}
\clearpage\appendix\section{Final estimates for every measured static quantity}
\input{tables/final-estimates.tex}
\clearpage\section{Representative plots for every retained case}
All physical static quantities, energy-moment audits and representative precision diagnostics have individual PDF, SVG and PNG assets in the figures folder. Raw-estimator audits are included in the tables and metadata. Exact references have zero sampling error; their numerical residuals are recorded in metadata.
\input{sections/figures.tex}
\end{document}
''',encoding='utf-8')
    (folder/'README.md').write_text(f'Open main.tex as the Overleaf main document. This is an organized folder, not a ZIP.\nStudy: {job["id"]}\nSource revision: {job["provenance"]["revision"]}\nPlots: {panels} panels in PDF/SVG/PNG. Exact data and failed cases: tables/ and metadata/.\nOriginal immutable sources, engines, inputs and traces: {STORE/job["id"]}\nLaTeX source/asset checks do not substitute for an actual Overleaf compilation.\nCopyright (c) 2026 Vikas Vijigiri. All rights reserved.\n',encoding='utf-8')
    manifest={p.relative_to(folder).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(folder.rglob('*')) if p.is_file() and p.name!='manifest.json'}
    (folder/'metadata/manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps(dict(folder=str(folder),panels=panels,hashed_files=len(manifest)),indent=2))


if __name__=='__main__':main()
