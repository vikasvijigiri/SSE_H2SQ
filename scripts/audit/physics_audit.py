"""Reproducible audit of geometry, exact classical traces, and QMC mixing."""
from pathlib import Path
import importlib.util
import sys, json, subprocess, itertools, time
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
from build_probes import prepare_baseline
prepare_baseline()
_spec = importlib.util.spec_from_file_location("historical_ed", ROOT / '.cache/physics-audit/historical-baseline/Analysis/ED_analysis/exact_diagonalzation_H2SQ.py')
ed = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = ed
_spec.loader.exec_module(ed)
OUT=ROOT/'.cache/physics-audit'
PROBE=OUT/'sse_probe.exe'

def geometry(lx,ly,lz=1):
    return json.loads(subprocess.check_output([str(PROBE),'geometry',str(lx),str(ly),str(lz)]))

def ed_terms(lx,ly):
    pl=ed.generate_black_plaquette_indices(lx,ly); nbr=ed.generate_plaquette_neighbor_map(lx,ly)
    return pl,nbr

def state_arrays(lx,ly):
    n=lx*ly
    bits=np.arange(1<<n,dtype=np.uint32)
    spins=2*((bits[:,None]>>np.arange(n))&1).astype(float)-1
    return bits,spins

def energies(lx,ly,j0,j1,j2,intended=False):
    bits,s=state_arrays(lx,ly); pl,nbr=ed_terms(lx,ly)
    energy=np.zeros(len(s)); qenergy=np.zeros(len(s)); g=geometry(lx,ly)
    for i,p in enumerate(pl):
        a,b,c,d=p.indices; aa=s[:,a]*s[:,c];bb=s[:,b]*s[:,d]
        energy-=j0*aa*bb
        energy+=j1*(aa+bb) if intended else abs(j1)*(1-aa-bb-aa*bb)
        for neighbor in nbr[i][:2]:
            q=pl[neighbor].indices
            px=(s[:,a]+s[:,b]-s[:,c]-s[:,d])/4
            py=(s[:,a]-s[:,b]+s[:,c]-s[:,d])/4
            qx=(s[:,q[0]]+s[:,q[1]]-s[:,q[2]]-s[:,q[3]])/4
            qy=(s[:,q[0]]-s[:,q[1]]+s[:,q[2]]-s[:,q[3]])/4
            energy+=j2*(px*qx+py*qy)
    # QMC site ordering x+Lx*y; transform to ED's x*Ly+y.
    mapped=np.array([(i%lx)*ly+i//lx for i in range(lx*ly)])
    qs=s[:,mapped]
    for p in g['plaquettes']:
        aa=qs[:,p[0]]*qs[:,p[2]];bb=qs[:,p[1]]*qs[:,p[3]]
        qenergy-=abs(j0)*aa*bb
        qenergy+=abs(j1)*(1-aa-bb-aa*bb)
    for a,b,c in g['j2bonds']:qenergy+=j2*c*qs[:,a]*qs[:,b]
    qenergy+=j2*g['j2constant']
    # Use the same physical checkerboard stagger under the ED indexing.
    phase=np.array([(-1)**(x+y) for x in range(lx) for y in range(ly)])
    mz=s@phase/(2*lx*ly)
    return energy,qenergy,mz

def exact(energy,mz,beta,n):
    w=np.exp(-beta*(energy-energy.min()));w/=w.sum(); e=energy/n
    return dict(E=float(w@e),E2=float(w@(e**2)),E4=float(w@(e**4)),m2=float(w@(mz**2)),m4=float(w@(mz**4)),Mx=0.0,ground_energy=float(energy.min()/n),ground_degeneracy=int(np.sum(np.isclose(energy,energy.min()))))

def run_case(beta,j1,seed,initial='hot',hx=0,j2=0,therm=50000,measure=105000):
    name=f'{PROBE.stem}_b{beta}_j{j1}_h{hx}_j2{j2}_seed{seed}_{initial}'
    command=[str(PROBE),str(beta),'1',str(j1),str(j2),str(hx),str(therm),str(measure),str(seed),initial,str(OUT/(name+'.csv'))]
    start=time.monotonic();result=json.loads(subprocess.check_output(command));result.update(beta=beta,J1=j1,J2=j2,hx=hx,seed=seed,initial=initial,seconds=time.monotonic()-start)
    print(json.dumps(result),flush=True)
    return result

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    report={'geometry':[],'exact':{},'runs':[]}
    for lx,ly in [(4,4),(4,2),(2,4)]:
        e,q,m=energies(lx,ly,1,.3,.025)
        report['geometry'].append(dict(lattice=[lx,ly,1],max_energy_difference=float(np.max(np.abs(e-q))),plaquettes=len(geometry(lx,ly)['plaquettes'])))
    for beta in [1,10]:
        for j1 in [0,.3]:
            e,q,m=energies(4,4,1,j1,0)
            report['exact'][f'code_beta{beta}_j1{j1}']=exact(e,m,beta,16)
            e,q,m=energies(4,4,1,j1,0,True)
            report['exact'][f'intended_beta{beta}_j1{j1}']=exact(e,m,beta,16)
    print(json.dumps({k:v for k,v in report.items() if k!='runs'},indent=2),flush=True)
    for beta,j1,seeds,initials in [(10,.3,[1,2,3],['hot']),(10,.3,[1],['cold','stripe']),(1,.3,[1,2,3],['hot']),(10,0,[1,2],['hot'])]:
        for initial in initials:
            for seed in seeds:
                report['runs'].append(run_case(beta,j1,seed,initial))
                (OUT/'audit.json').write_text(json.dumps(report,indent=2))
    (OUT/'audit.json').write_text(json.dumps(report,indent=2))

if __name__=='__main__':main()
