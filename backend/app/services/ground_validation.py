# Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
"""Reproducible random-coupling ground-state validation with explicit NO GO cases."""
import asyncio
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import hashlib
import io
import json
import math
import platform
from pathlib import Path
import subprocess
import threading
import time
import uuid

import numpy as np
import pandas as pd
from pydantic import BaseModel,ConfigDict,Field,model_validator
from scipy.stats import poisson
from . import ground_ed
from .ground_statistics import summarize
from .physics_model import exact_reference,MODEL_VERSION

ROOT=Path(__file__).resolve().parents[3]
STORE=ROOT/'.cache/ground-validation'
LOCK=threading.RLock();ACTIVE={}
PHYSICS_KEYS=['E','mz2','mz4','mx','plaquette','ice','rho','stripe2','stripe4','polar2','polar4','ct0']
PROTOCOL='open-projector-random-validation-v1'
PYTHON_SOURCES=[Path(__file__).with_name(name) for name in ['ground_ed.py','ground_statistics.py','ground_validation.py','physics_model.py','benchmark_statistics.py','symmetry_ed.py']]
LOADED_PYTHON_SHA={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in PYTHON_SOURCES}


class Request(BaseModel):
    model_config=ConfigDict(allow_inf_nan=False)
    Lx:int=4;Ly:int=4
    points:int=Field(default=8,ge=1,le=24)
    seed:int=Field(default=20261021,ge=0,le=2000000000)
    chains:int=Field(default=8,ge=4,le=16)
    n_therm:int=Field(default=30000,ge=1000,le=1000000)
    n_measure:int=Field(default=105000,ge=10000,le=1000000)
    lower_measure:int=Field(default=35000,ge=10000,le=1000000)
    max_projection:float=Field(default=256,ge=8,le=4096)
    fidelity_tolerance:float=Field(default=1e-10,ge=1e-12,le=1e-4)
    relative_precision:float=Field(default=.01,gt=0,le=.1)
    histogram_fidelity_min:float=Field(default=.995,gt=.9,lt=1.)
    histogram_tvd_max:float=Field(default=.025,gt=0,le=.2)
    histogram_missing_weight_max:float=Field(default=.005,gt=0,le=.1)
    include_controls:bool=True
    parameters:list[dict[str,float]]|None=None

    @model_validator(mode='after')
    def valid(self):
        if min(self.Lx,self.Ly)<2 or self.Lx%2 or self.Ly%2 or self.Lx*self.Ly>256:
            raise ValueError('Use even periodic 2D cells; this machine profile permits at most 256 spins.')
        if self.parameters is not None:
            if not 1<=len(self.parameters)<=24:raise ValueError('Provide 1-24 parameter sets.')
            for p in self.parameters:
                if set(p)!={'J0','J1','J2','hx'} or not all(np.isfinite(v) and abs(v)<=100 for v in p.values()) or p['hx']==0:
                    raise ValueError('Each quantum set needs finite J0,J1,J2,hx and nonzero hx. h=0 is covered by the separate classical check.')
        return self


def save(job):
    with LOCK:
        folder=STORE/job['id'];folder.mkdir(parents=True,exist_ok=True)
        temporary=folder/'study.tmp';temporary.write_text(json.dumps(job,allow_nan=False),encoding='utf-8')
        for attempt in range(20):
            try:temporary.replace(folder/'study.json');break
            except PermissionError:
                if attempt==19:raise
                time.sleep(.05)


def get(identity):
    if len(identity)!=32 or any(c not in '0123456789abcdef' for c in identity):raise FileNotFoundError(identity)
    with LOCK:
        job=json.loads((STORE/identity/'study.json').read_text(encoding='utf-8'))
        if job['status'] in ('running','queued') and identity not in ACTIVE:
            job.update(status='interrupted',message='Service restarted; original traces remain available. Resume after source verification.');save(job)
        return job


def list_jobs():
    return sorted([{k:get(p.parent.name)[k] for k in ('id','created','status','completed_runs','total_runs','verdict')} for p in STORE.glob('*/study.json')],key=lambda j:j['created'],reverse=True)[:20]


def provenance():
    paths=[ROOT/'backend/build/Makefile',ROOT/'backend/benchmarks/projector.cpp',ROOT/'backend/benchmarks/emc.cpp',ROOT/'backend/benchmarks/runner.cpp',*sorted((ROOT/'backend/src').glob('*.cpp')),*sorted((ROOT/'backend/include').glob('*.hpp')),*[Path(__file__).with_name(name) for name in ['ground_ed.py','ground_statistics.py','ground_validation.py','physics_model.py','benchmark_statistics.py','symmetry_ed.py']]]
    return dict(protocol=PROTOCOL,model_version=MODEL_VERSION,revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True,cwd=ROOT).strip(),working_tree_dirty=bool(subprocess.check_output(['git','status','--porcelain'],cwd=ROOT)),
        source_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        engine_sha256={name:hashlib.sha256((ROOT/'backend/build'/f'{name}.exe').read_bytes()).hexdigest() for name in ['sse_projector','emc_validation','sse_benchmark']},
        compiler=subprocess.check_output(['g++','--version'],text=True).splitlines()[0],flags='-O2 -std=c++17 -Wall; no fast-math',python=platform.python_version(),numpy=np.__version__)


def sample_parameters(request):
    if request.parameters is not None:return request.parameters
    rng=np.random.default_rng(request.seed);order=rng.permutation(8);result=[]
    for i in range(request.points):
        octant=int(order[i%8]);magnitudes=rng.uniform([.12,.025,.1],[.8,.15,.8])
        result.append(dict(J0=1.,J1=float(magnitudes[0]*(1 if octant&1 else -1)),J2=float(magnitudes[1]*(1 if octant&2 else -1)),hx=float(magnitudes[2]*(1 if octant&4 else -1))))
    return result


def analytic_ground(lx,ly,h):
    n=lx*ly;values=dict(E=-abs(h),E_count=-abs(h),E2=h*h,E4=h**4,mz2=1/(4*n),mz4=(3*n-2)/(16*n**3),mx=float(np.sign(h)),plaquette=0.,ice=.375,rho=1/6,stripe2=1/(2*n),stripe4=(n-1)/(2*n**3),polar2=2/n,polar4=8*(n-1)/n**3,ct0=.25)
    values.update({key:values[original] for key,original in ground_ed.RAW_ALIASES.items()})
    return dict(values=values,ground_energy_per_spin=-abs(h),ground_fidelity=1.,ground_infidelity=0.,numerical_residual=0.,dimension=2**n,block_dimension=None,projection_total=None,trial='exact product |sign(h)x>',reference='analytic noninteracting transverse-field ground state',sampling_error=0.,sector_gap=2*abs(h),energy_variance=0.)


def plan(job,request):
    points=[];base=dict(Lx=request.Lx,Ly=request.Ly)
    configs=[dict(label=f'{"Custom" if request.parameters is not None else "Random"} {i+1}',kind='ground',params={**base,**p}) for i,p in enumerate(job['sampled_parameters'])]
    if request.include_controls:
        configs+=[dict(label='Screenshot ground-state challenge',kind='ground',params=dict(Lx=4,Ly=4,J0=1.,J1=.3,J2=.025,hx=.2)),
                  dict(label='Near-zero coupling stress; no automatic convergence claim',kind='stress',params=dict(Lx=4,Ly=4,J0=1.,J1=.3,J2=.0001,hx=.02))]
        for i,(lx,ly) in enumerate([(2,2),(4,2),(2,6),(2,8),(8,2)]):configs.append(dict(label=f'{lx}x{ly} independent-size probe',kind='ground',params=dict(Lx=lx,Ly=ly,**job['sampled_parameters'][i%len(job['sampled_parameters'])])))
        for lx,ly in [(4,8),(6,6),(8,8)]:configs.append(dict(label=f'{lx}x{ly} exact product-ground limit',kind='analytic',params=dict(Lx=lx,Ly=ly,J0=0.,J1=0.,J2=0.,hx=.37)))
    for index,item in enumerate(configs):
        job['message']=f'Planning exact reference: {item["label"]}';save(job)
        p=item['params'];item.update(id=f'p{index:03d}',levels=[],verdict='PENDING')
        try:
            if item['kind']=='analytic':
                item['ground']=analytic_ground(p['Lx'],p['Ly'],p['hx']);lengths=[2.,4.,8.]
                for length in lengths:item['levels'].append(dict(projection=length,ed={**item['ground'],'projection_total':length},methods={}))
            else:
                ground=ground_ed.reference(p);item['ground']=ground
                if ground['sector_gap']<100*ground['numerical_residual']:
                    raise RuntimeError('Ground sector gap is not numerically resolved relative to the ED residual.')
                if math.sqrt(ground['block_dimension'])*ground['numerical_residual']/ground['sector_gap']>math.sqrt(request.fidelity_tolerance)/10:
                    raise RuntimeError('Numerical ground-eigenvector resolution is inadequate for the requested fidelity target.')
                grid=[2.**k for k in range(2,int(math.floor(math.log2(request.max_projection)))+1)]
                refs=[ground_ed.reference(p,t) for t in grid]
                eligible=[i for i,r in enumerate(refs) if r['ground_infidelity']<=request.fidelity_tolerance]
                terminal=min(eligible[0]+1,len(grid)-1) if eligible else len(grid)-1
                lengths=[grid[terminal]/4,grid[terminal]/2,grid[terminal]] if item['kind']!='stress' else [grid[-1]]
                for length in lengths:item['levels'].append(dict(projection=length,ed=ground_ed.reference(p,length),methods={}))
                item['projection_plan_converged']=bool(eligible)
                if not eligible:item['planning_warning']='Ground fidelity target not reached within the projection budget; NO GO for a ground-state claim.'
        except (ValueError,RuntimeError) as error:
            item.update(verdict='NO GO',reference_error=str(error))
        if item.get('ground') and item['kind'] in ('ground','stress') and item['levels']:
            directory=STORE/job['id']/'references'/item['id'];directory.mkdir(parents=True,exist_ok=True)
            mapping,probability=ground_ed.orbit_probabilities(p)
            map_path=directory/'orbit-map.bin';mapping.tofile(map_path)
            ground_path=directory/'ground-probability.npy';np.save(ground_path,probability)
            item['orbit_reference']=dict(map_file=map_path.relative_to(STORE/job['id']).as_posix(),ground_file=ground_path.relative_to(STORE/job['id']).as_posix(),
                orbits=len(probability),map_sha256=hashlib.sha256(map_path.read_bytes()).hexdigest(),ground_sha256=hashlib.sha256(ground_path.read_bytes()).hexdigest())
            for level in item['levels']:
                _,probability=ground_ed.orbit_probabilities(p,level['projection']);path=directory/f'projected-{level["projection"]:g}.npy';np.save(path,probability)
                level['orbit_probability_file']=path.relative_to(STORE/job['id']).as_posix();level['orbit_probability_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
        points.append(item);job['points']=points;save(job)
    if request.include_controls:
        for index,beta in enumerate([.5,2.,8.]):
            p=dict(Lx=4,Ly=4,**job['sampled_parameters'][index%len(job['sampled_parameters'])]);p['hx']=0.;p['beta']=beta
            ed=exact_reference(p);ed['values']['E_count']=ed['values']['E'];ed['values'].update({key:ed['values'][original] for key,original in ground_ed.RAW_ALIASES.items()});ed['values']={k:ed['values'][k] for k in ground_ed.KEYS};ed['ground_fidelity']=None
            points.append(dict(id=f'c{index:03d}',label=f'Independent classical EMC / SSE / exact trace, beta={beta:g}',kind='classical',params=p,ground=ed,levels=[dict(projection=beta,ed=ed,methods={})],verdict='PENDING'))
        lx,ly=(request.Lx,request.Ly) if request.Lx*request.Ly>16 else (6,6)
        p=dict(Lx=lx,Ly=ly,**job['sampled_parameters'][0]);p.update(hx=0.,beta=.5)
        unknown=dict(values={key:0. if key in ('mx','raw_mx') else None for key in ground_ed.KEYS},ground_fidelity=None,numerical_residual=None,
                     reference='independent EMC/SSE cross-check; exact interacting reference unavailable',dimension=2**(p['Lx']*p['Ly']))
        points.append(dict(id='large-classical',label=f'{p["Lx"]}x{p["Ly"]} independent classical cross-check; no ED claim',kind='classical_cross',params=p,ground=unknown,levels=[dict(projection=.5,ed=unknown,methods={})],verdict='PENDING'))
    for index,p in enumerate(points):p['seed_offset']=index*100000
    job['points']=points
    job['total_runs']=sum(len(p['levels'])*request.chains*(2 if p['kind'].startswith('classical') else 1) for p in points)
    job['planning_complete']=True;save(job)
    (STORE/job['id']/'plan.json').write_text(json.dumps(dict(request=job['request'],sampled_parameters=job['sampled_parameters'],points=points,provenance=job['provenance']),allow_nan=False),encoding='utf-8')


def run_chain(job,point,level,method,chain,request,cancelled):
    if cancelled.is_set():raise InterruptedError('Cancelled; completed original traces retained.')
    final=level is point['levels'][-1]
    measurements=request.n_measure if final and point['kind']!='stress' else request.lower_measure
    params={**point['params'],'beta':level['projection'],'Lz':1,'J3':0,'n_therm':request.n_therm,'n_measure':measurements}
    directory=STORE/job['id']/point['id']/f'{level["projection"]:g}'/method/str(chain);directory.mkdir(parents=True,exist_ok=True)
    text='\n'.join(f'{k} {v}' for k,v in params.items())+'\n'
    seed=request.seed+point['seed_offset']+point['levels'].index(level)*10000+{'projector':0,'emc':1000,'sse':2000}[method]+chain
    if (directory/'trace.npz').exists() and (directory/'metadata.json').exists() and (directory/'input.in').read_text()==text:
        with np.load(directory/'trace.npz',allow_pickle=False) as saved:trace=saved['trace']
        metadata=json.loads((directory/'metadata.json').read_text())
        if len(trace)==measurements and metadata['seed']==seed:return trace,metadata
    (directory/'input.in').write_text(text)
    executable={'projector':'sse_projector','emc':'emc_validation','sse':'sse_benchmark'}[method]
    command=[str(STORE/job['id']/'engines'/f'{executable}.exe'),str(directory/'input.in'),str(directory/'trace.csv')]
    if hashlib.sha256(Path(command[0]).read_bytes()).hexdigest()!=job['provenance']['engine_sha256'][executable]:
        raise RuntimeError('Frozen validation executable changed.')
    command += ['symmetry',str(seed)] if method=='sse' else [str(seed),'run'] if method=='emc' else [str(seed)]
    if method=='projector' and point.get('orbit_reference'):
        map_path=STORE/job['id']/point['orbit_reference']['map_file']
        if hashlib.sha256(map_path.read_bytes()).hexdigest()!=point['orbit_reference']['map_sha256']:raise RuntimeError('Frozen orbit map changed.')
        command.append(str(map_path))
    started=time.monotonic();process=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,cwd=ROOT)
    try:
        while True:
            try:stdout,stderr=process.communicate(timeout=.5);break
            except subprocess.TimeoutExpired:
                if cancelled.is_set():process.terminate();process.communicate();raise InterruptedError('Cancelled; completed original traces retained.')
        if process.returncode:raise RuntimeError(f'{method} failed: {stderr[-1000:]}')
    finally:
        if process.poll() is None:process.terminate();process.communicate()
    metadata=json.loads(stdout);metadata.update(seed=seed,total_wall_seconds=time.monotonic()-started,measurement_rows=measurements)
    data=pd.read_csv(directory/'trace.csv');trace=np.empty(len(data),dtype=[(key,'f8') for key in ground_ed.KEYS])
    for key in ground_ed.KEYS:trace[key]=data[key if key in data else 'E' if key=='E_count' else ground_ed.RAW_ALIASES[key]]
    if not all(np.isfinite(trace[key]).all() for key in trace.dtype.names) or len(trace)!=measurements:raise RuntimeError('Invalid trace rows.')
    arrays=dict(trace=trace)
    if method=='projector' and point.get('orbit_reference'):
        histogram=np.fromfile(directory/'trace.csv.hist.bin',dtype='<f8')
        if len(histogram)!=point['orbit_reference']['orbits'] or np.any(histogram<0) or not np.isclose(histogram.sum(),len(trace),atol=1e-6,rtol=0):
            raise RuntimeError('Invalid midpoint orbit-histogram normalization.')
        arrays['histogram']=histogram;(directory/'trace.csv.hist.bin').unlink()
    np.savez_compressed(directory/'trace.npz',**arrays);(directory/'metadata.json').write_text(json.dumps(metadata));(directory/'trace.csv').unlink()
    return trace,metadata


def evaluate_point(point,request):
    if not point['levels']:return
    final=point['levels'][-1];method='sse' if point['kind'].startswith('classical') else 'projector'
    summary=final['methods'][method];obs=summary['observables'];reasons=[]
    if any(v['flags'] for v in obs.values()):reasons.append('sampling or finite-ensemble ED checks failed')
    if any(not obs[key]['precision_pass'] for key in PHYSICS_KEYS):reasons.append('requested statistical precision not reached')
    if point['kind'].startswith('classical'):
        other=final['methods']['emc'];differences={}
        for key,v in obs.items():
            w=other['observables'][key];se=math.hypot(v['error'] or 0,w['error'] or 0);delta=v['mean']-w['mean'];differences[key]=dict(difference=delta,combined_SE=se,pull=delta/se if se else None)
            if abs(delta)>4*se+1e-10:reasons.append(f'EMC/SSE disagreement: {key}')
        if any(v['flags'] for v in other['observables'].values()):reasons.append('independent EMC convergence checks failed')
        point['independent_comparison']=differences
        point['verdict_scope']='finite-temperature classical equilibrium; not quantum ground-state validation'
    else:
        fidelity=final['ed']['ground_fidelity']
        if final['ed'].get('ground_infidelity',1-fidelity)>request.fidelity_tolerance:reasons.append('pure-state projection fidelity target not reached')
        if any(not obs[key]['projection_bias_pass'] or not obs[key]['ground_agreement'] for key in PHYSICS_KEYS):reasons.append('ground-state observable convergence failed')
        for level in point['levels']:
            if any(v['flags'] for v in level['methods']['projector']['observables'].values()):reasons.append(f'finite-projector check failed at T={level["projection"]:g}')
        if point.get('orbit_reference'):
            histogram=summary.get('histogram')
            if not histogram:reasons.append('full orbit-distribution comparison missing')
            elif histogram['fidelity_ground']<request.histogram_fidelity_min or histogram['tvd_ground']>request.histogram_tvd_max or histogram['missing_ground_weight']>request.histogram_missing_weight_max:
                reasons.append('sampled ground-state distribution tolerance not reached')
            elif any(r['fidelity_projected']<.95 or r['tvd_projected']>.2 for r in histogram['per_run']):
                reasons.append('independent-run orbit populations require more equilibration or statistics')
        if len(point['levels'])>=2:
            if point['levels'][-2]['ed'].get('ground_infidelity',1-point['levels'][-2]['ed']['ground_fidelity'])>request.fidelity_tolerance:
                reasons.append('two final projection lengths have not both reached the state-fidelity target')
            previous=point['levels'][-2]['methods']['projector']['observables'];plateau={}
            for key in PHYSICS_KEYS:
                error=math.hypot(obs[key]['error'] or 0,previous[key]['error'] or 0)
                delta=obs[key]['mean']-previous[key]['mean'];exact_delta=final['ed']['values'][key]-point['levels'][-2]['ed']['values'][key]
                passed=abs(delta-exact_delta)<=4*error+1e-10 and abs(delta)<=4*error+obs[key]['precision_tolerance']/4
                plateau[key]=dict(difference=delta,combined_SE=error,exact_projection_difference=exact_delta,passed=passed)
                if not passed:reasons.append(f'projection-length plateau unresolved: {key}')
            point['plateau']=plateau
        else:reasons.append('projection-length plateau not measured')
        point['verdict_scope']='pure-state static ground observables at the stated lattice, couplings, precision and projection budget'
    point['verdict']='NO GO' if reasons else 'GO';point['reasons']=sorted(set(reasons))


def cutoff_bias_bounds(runs,params,method):
    """Tail and count-moment controls also apply to positive projector overlaps."""
    bounds=[];n=params['Lx']*params['Ly'];h=abs(params['hx'])
    diagonal_norm=.5*abs(params['J0'])+abs(params['J1'])+abs(params['J2'])
    norms=dict(mz2=.25,mz4=.0625,plaquette=1.,ice=1.,rho=1.,stripe2=.5,stripe4=.25,polar2=2.,polar4=4.,ct0=1.)
    for r in runs:
        m,mu=r['operator_cutoff'],r['norm_mean_bound'];tail=float(poisson.sf(m,mu));total=r['inverse_temperature'];shift=abs(r['energy_shift_total'])
        if tail>=1:raise RuntimeError('SSE cutoff discards the whole probability bound.')
        delta={k:mu**k*(float(poisson.sf(m-k,mu))+tail)/(1-tail) for k in range(1,5)}
        values={key:2*norm*tail/(1-tail) for key,norm in norms.items()}
        for power,key in [(1,'E_count'),(2,'E2'),(4,'E4')]:
            values[key]=sum(math.comb(power,k)*shift**(power-k)*delta[k]/total**k for k in range(1,power+1))/n**power
        values['mx']=delta[1]/(h*total*n) if h else 0.
        values.update({key:values[original] for key,original in ground_ed.RAW_ALIASES.items()})
        values['E']=2*diagonal_norm*tail/(1-tail)+delta[1]/(total*n) if method=='projector' else values['E_count']
        r['sse_partition_tail_bound']=tail;r['sse_observable_bias_bounds']=values;bounds.append(values)
    return {key:max(v[key] for v in bounds) for key in ground_ed.KEYS}


def histogram_summary(job,point,level,method,runs):
    folder=STORE/job['id'];histograms=[]
    for chain in range(len(runs)):
        with np.load(folder/point['id']/f'{level["projection"]:g}'/method/str(chain)/'trace.npz',allow_pickle=False) as saved:histograms.append(saved['histogram'])
    for path,digest in [(folder/level['orbit_probability_file'],level['orbit_probability_sha256']),(folder/point['orbit_reference']['ground_file'],point['orbit_reference']['ground_sha256'])]:
        if hashlib.sha256(path.read_bytes()).hexdigest()!=digest:raise RuntimeError('Frozen reference probabilities changed.')
    histograms=np.asarray(histograms);q=np.load(folder/level['orbit_probability_file'],allow_pickle=False);g=np.load(folder/point['orbit_reference']['ground_file'],allow_pickle=False)
    def metrics(counts):
        p=counts/counts.sum()
        return dict(fidelity_projected=float(np.clip(np.sum(np.sqrt(p*q))**2,0,1)),fidelity_ground=float(np.clip(np.sum(np.sqrt(p*g))**2,0,1)),
            tvd_projected=float(np.abs(p-q).sum()/2),tvd_ground=float(np.abs(p-g).sum()/2),
            missing_projected_weight=float(q[p==0].sum()),missing_ground_weight=float(g[p==0].sum()))
    counts=histograms.sum(0);result=metrics(counts);rng=np.random.default_rng(8271964)
    bootstrap=[metrics(histograms[rng.integers(0,len(histograms),len(histograms))].sum(0)) for _ in range(200)]
    result['bootstrap_SE']={key:float(np.std([b[key] for b in bootstrap],ddof=1)) for key in result}
    result['per_run']=[metrics(counts) for counts in histograms]
    result['prefix_history']=[dict(runs=count,rows=float(histograms[:count].sum()),**metrics(histograms[:count].sum(0))) for count in sorted(set([min(4,len(histograms)),len(histograms)]))]
    result['orbits']=len(q);result['rows']=float(counts.sum())
    result['reconstructed_state']=ground_ed.reconstructed_diagnostics(point['params'],counts)
    result['uncertainty_note']='Whole independent-run bootstrap (200 replicates), conditional on observed support. Missing exact-reference weight is reported; these SEs are descriptive, not a 100% coverage guarantee.'
    return result


def execute(job,request,cancelled):
    try:
        job.update(status='running',message='Freezing random parameter list and exact projection plan');save(job)
        if not job.get('planning_complete'):plan(job,request)
        job['completed_runs']=0;save(job)
        for point in job['points']:
            for level in point['levels']:
                for method in (['emc','sse'] if point['kind'].startswith('classical') else ['projector']):
                    job['message']=f'{point["label"]}; {method}; total projection / beta={level["projection"]:g}';save(job)
                    completed={}
                    with ThreadPoolExecutor(max_workers=min(4,request.chains)) as pool:
                        futures={pool.submit(run_chain,job,point,level,method,i,request,cancelled):i for i in range(request.chains)}
                        for future in as_completed(futures):
                            completed[futures[future]]=future.result();job['completed_runs']+=1;save(job)
                    traces,runs=zip(*(completed[i] for i in range(request.chains)))
                    seconds=sum(r['total_wall_seconds'] for r in runs)
                    constants=['mx','raw_mx'] if point['kind'].startswith('classical') else []
                    if method=='projector' and all(point['params'][key]==0 for key in ('J0','J1','J2')):
                        constants+=['mz2','mz4','plaquette','ice','rho','stripe2','stripe4','polar2','polar4','ct0']
                    summary=summarize(traces,level['ed'],point['ground'],seconds,request.relative_precision,constants)
                    summary['runs']=runs;summary['sampling_checks']=[]
                    summary['estimator_comparison']={original:dict(raw_SE=summary['observables'][key]['error'],conditional_SE=summary['observables'][original]['error'],
                        SE_gain=summary['observables'][key]['error']/summary['observables'][original]['error'] if summary['observables'][key]['error'] is not None and summary['observables'][original]['error'] else None) for key,original in ground_ed.RAW_ALIASES.items()}
                    summary['update_seconds']={stage:sum(r.get('update_seconds',{}).get(stage,0.) for r in runs) for stage in ['diagonal','cluster','gauge','canonical','symmetry']}
                    if method!='emc':
                        bounds=[float(poisson.sf(r['operator_cutoff'],r['norm_mean_bound'])) for r in runs]
                        summary['sse_partition_tail_bound']=max(bounds)
                        summary['sse_observable_bias_bounds']=cutoff_bias_bounds(runs,point['params'],method)
                        for key,v in summary['observables'].items():v['sse_absolute_bias_bound']=summary['sse_observable_bias_bounds'][key]
                        if max(bounds)>1e-20 or any(r['cutoff_saturations'] for r in runs):summary['sampling_checks'].append('SSE cutoff bound or saturation requires review')
                        if method=='projector' and any(r['open_endpoint_rows']==0 for r in runs):summary['sampling_checks'].append('open projector endpoints were never observed')
                        if method=='projector' and any(r.get('local_skip_probability',0)<=0 or len(r.get('site_symmetry_orbit_counts',[]))!=4 or sum(r['site_symmetry_orbit_counts'])!=r['measurement_rows'] or min(r['site_symmetry_orbit_counts'])==0 for r in runs):
                            summary['sampling_checks'].append('projector identity or symmetry-orbit activity missing')
                    elif any(r['round_trips']<10 for r in runs):summary['sampling_checks'].append('fewer than ten measured temperature round trips in EMC')
                    if summary['sampling_checks']:
                        for v in summary['observables'].values():v['flags']+=summary['sampling_checks'];v['status']='NO GO'
                    if len(traces)>=8 and level is point['levels'][-1]:
                        first=summarize(traces[:4],level['ed'],point['ground'],sum(r['total_wall_seconds'] for r in runs[:4]),request.relative_precision,constants)
                        history=[]
                        for count,estimates in [(4,first),(len(traces),summary)]:
                            history.append(dict(runs=count,rows=sum(len(t) for t in traces[:count]),expected_SE_ratio=math.sqrt(sum(len(t) for t in traces[:4])/sum(len(t) for t in traces[:count])),
                                observables={key:{**{name:v[name] for name in ['mean','error','tau_int','ess','rhat','status']},'SE_ratio':v['error']/first['observables'][key]['error'] if v['error'] is not None and first['observables'][key]['error'] else None} for key,v in estimates['observables'].items()}))
                        summary['precision_history']=history
                    if method=='projector' and point.get('orbit_reference'):
                        summary['histogram']=histogram_summary(job,point,level,method,runs)
                    level['methods'][method]=summary;save(job)
            evaluate_point(point,request);save(job)
        random_points=[p for p in job['points'] if p['label'].startswith(('Random','Custom'))]
        job['verdict'].update(random_ground_domain='GO' if random_points and all(p['verdict']=='GO' for p in random_points) else 'NO GO',
            completed_ground_points=sum(p['verdict']=='GO' and p['kind'] in ('ground','analytic') for p in job['points']),
            failed_points=[p['id'] for p in job['points'] if p['verdict']!='GO'])
        job.update(status='complete',message='Validation finished. Read per-point verdicts; the universal all-lattice/all-coupling claim remains NO GO.')
    except InterruptedError as error:job.update(status='cancelled',message=str(error))
    except Exception as error:job.update(status='failed',message=str(error))
    finally:
        job['finished']=datetime.now(timezone.utc).isoformat();save(job)
        with LOCK:ACTIVE.pop(job['id'],None)


async def start(request):
    from .benchmark import ACTIVE as thermal_active
    if any(hashlib.sha256(p.read_bytes()).hexdigest()!=LOADED_PYTHON_SHA[p.name] for p in PYTHON_SOURCES):
        raise RuntimeError('Python numerical source changed since service startup. Run DEPLOY ALL before freezing a study.')
    dependencies=[*sorted((ROOT/'backend/include').glob('*.hpp')),ROOT/'backend/src/updates.cpp',ROOT/'backend/src/variables.cpp',ROOT/'backend/src/observables.cpp']
    for engine,source in [('sse_projector','projector.cpp'),('emc_validation','emc.cpp'),('sse_benchmark','runner.cpp')]:
        executable=ROOT/'backend/build'/f'{engine}.exe'
        if not executable.exists() or any(p.stat().st_mtime>executable.stat().st_mtime for p in [*dependencies,ROOT/'backend/benchmarks'/source]):
            raise RuntimeError('Validation engine is missing or stale. Run DEPLOY ALL first.')
    with LOCK:
        if ACTIVE or thermal_active:raise RuntimeError('Another validation study is active.')
        identity=uuid.uuid4().hex;event=threading.Event();ACTIVE[identity]=event
    try:
        info=provenance();job=dict(id=identity,created=datetime.now(timezone.utc).isoformat(),status='queued',request=request.model_dump(),
            sampled_parameters=sample_parameters(request),points=[],completed_runs=0,total_runs=0,provenance=info,
            verdict=dict(universal='NO GO',random_ground_domain='PENDING',reason='A finite dataset cannot certify every real coupling and every lattice size.'),
            definitions=dict(ensemble='|psi(T/2)> proportional to exp(-T H/2)|+x>; open imaginary-time ends; no thermal trace',
                quantities='Static diagonal centre measurements, exact conditional off-diagonal event density at centre, count-energy audits; no zero-temperature susceptibility claim',
                scope='Even periodic 2D checkerboard cells. Exact interacting references are resource-limited; unsupported cells remain NO GO.',
                randomness='User-supplied coupling sets; not random samples.' if request.parameters is not None else 'PCG64 seeded in advance; signed octants balanced; |J1| in [0.12,0.8], |J2| in [0.025,0.15], |h| in [0.1,0.8]; all couplings vary simultaneously.',
                distributions='On ED-sized quantum cells, exact symmetry-orbit averages reconstruct the full positive-state z distribution. Raw Bhattacharyya fidelity, total variation, missing ED weight, independent-run bootstrap and variational energy are reported. No ED-informed pseudocounts are added.',
                go='Finite-ensemble ED checks, ground fidelity, projection-length plateau, correlated precision, R-hat, ESS, block adequacy and cutoff checks must all pass.',
                confidence='One-SE bars and four-SE review criteria are statistical diagnostics, not a 100% coverage or universal mixing guarantee.'))
        save(job)
        for name,digest in info['source_sha256'].items():
            data=(ROOT/name).read_bytes();assert hashlib.sha256(data).hexdigest()==digest
            path=STORE/identity/'sources'/digest/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
        for name,digest in info['engine_sha256'].items():
            data=(ROOT/'backend/build'/f'{name}.exe').read_bytes();assert hashlib.sha256(data).hexdigest()==digest
            path=STORE/identity/'engines'/f'{name}.exe';path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
        asyncio.create_task(asyncio.to_thread(execute,job,request,event));return job
    except Exception:
        with LOCK:ACTIVE.pop(identity,None)
        raise


def cancel(identity):
    with LOCK:
        if identity in ACTIVE:ACTIVE[identity].set();return True
    return False


async def resume(identity):
    job=get(identity)
    with LOCK:
        if ACTIVE:raise RuntimeError('Another ground validation is active.')
        if job['status'] not in ('cancelled','failed','interrupted'):raise RuntimeError('Only stopped studies can resume.')
        current=provenance()
        if current['source_sha256']!=job['provenance']['source_sha256'] or current['engine_sha256']!=job['provenance']['engine_sha256']:
            raise RuntimeError('Numerical source or engine changed. Preserve this study and start a new frozen cohort.')
        event=threading.Event();ACTIVE[identity]=event
    job.update(status='queued',message='Resuming verified original traces');save(job)
    asyncio.create_task(asyncio.to_thread(execute,job,Request(**job['request']),event));return job


def trace_path(job,point_id,projection,method,chain):
    point=next((p for p in job['points'] if p['id']==point_id),None)
    if point is None or not any(l['projection']==projection and method in l['methods'] for l in point['levels']) or not 0<=chain<job['request']['chains']:
        raise FileNotFoundError('Unknown retained trace')
    return STORE/job['id']/point_id/f'{projection:g}'/method/str(chain)/'trace.npz'


def csv_result(job):
    import csv
    output=io.StringIO();writer=csv.writer(output);writer.writerow(['point','kind','Lx','Ly','J0','J1','J2','hx','projection_total_or_beta','method','observable','QMC_mean','SE','projected_ED','ground_ED','fidelity','infidelity','projection_bias','SSE_absolute_bias_bound','precision_tolerance','precision_pass','tau_int','ESS','Rhat','point_verdict','flags'])
    for p in job['points']:
        for l in p['levels']:
            for method,summary in l['methods'].items():
                for key,v in summary['observables'].items():writer.writerow([p['id'],p['kind'],*[p['params'][k] for k in ['Lx','Ly','J0','J1','J2','hx']],l['projection'],method,key,v['mean'],v['error'],v['projected_ed'],v['ground_ed'],l['ed'].get('ground_fidelity'),l['ed'].get('ground_infidelity'),v['projection_bias'],v.get('sse_absolute_bias_bound'),v['precision_tolerance'],v['precision_pass'],v['tau_int'],v['ess'],v['rhat'],p['verdict'],'; '.join(v['flags'])])
    return output.getvalue()


def figure_result(job,point_id,observable,quantity='mean',format='pdf'):
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    import matplotlib as mpl
    mpl.rcParams.update({'pdf.fonttype':42,'svg.fonttype':'none','font.family':'serif','mathtext.fontset':'stix'})
    point=next(p for p in job['points'] if p['id']==point_id)
    if quantity=='fidelity' and all(l['ed'].get('ground_fidelity') is None for l in point['levels']):
        raise ValueError('This classical cross-check has no quantum-state fidelity reference.')
    if quantity.startswith('state_') and not any(l['methods'].get('projector',{}).get('histogram') for l in point['levels']):
        raise ValueError('Configuration-distribution measurements are not available for this point.')
    figure=Figure(figsize=(6.4,4.4),layout='constrained');FigureCanvasAgg(figure)
    figure.get_layout_engine().set(rect=(0,.07,1,.92));ax=figure.subplots();ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.15)
    levels=point['levels'];x=[l['projection'] for l in levels]
    exact=[l['ed']['values'].get(observable) for l in levels];ground=point.get('ground',{}).get('values',{}).get(observable)
    if quantity=='mean':
        if any(v is not None for v in exact):ax.plot(x,exact,'D--',color='#6b7280',markersize=4,label='Exact finite ensemble')
        if ground is not None:ax.axhline(ground,color='black',linestyle=':',label='Ground ED / analytic' if not point['kind'].startswith('classical') else 'Exact classical trace')
    elif quantity in ('difference','sampling_difference'):ax.axhline(0,color='black',linestyle='--')
    if quantity.startswith('state_'):
        prefix={'state_infidelity':'fidelity','state_tvd':'tvd','state_missing':'missing'}[quantity]
        for reference,color,marker in [('projected','#0072B2','o'),('ground','#D55E00','s')]:
            key=f'{prefix}_{reference}'+('_weight' if prefix=='missing' else '')
            rows=[(l['projection'],l['methods']['projector']['histogram']) for l in levels if l['methods'].get('projector',{}).get('histogram')]
            y=[1-h[key] if prefix=='fidelity' else h[key] for _,h in rows]
            ax.errorbar([v for v,_ in rows],y,yerr=[h['bootstrap_SE'][key] for _,h in rows],fmt=marker+'-',color=color,capsize=3,label=f'QMC distribution vs {reference} ED')
        target={'state_infidelity':1-job['request']['histogram_fidelity_min'],'state_tvd':job['request']['histogram_tvd_max'],'state_missing':job['request']['histogram_missing_weight_max']}[quantity]
        ax.axhline(target,color='black',linestyle='--',label='Final ground review threshold')
    elif quantity=='fidelity':
        y=[max(1e-16,l['ed'].get('ground_infidelity',1-l['ed']['ground_fidelity'])) if l['ed'].get('ground_fidelity') is not None else None for l in levels]
        ax.plot(x,y,'D-',color='black',label='Exact projector infidelity');ax.axhline(job['request']['fidelity_tolerance'],color='#D55E00',linestyle='--',label='Required target');ax.set_yscale('log')
    else:
        for method,color,marker in [('projector','#0072B2','o'),('emc','#D55E00','s'),('sse','#009E73','^')]:
            rows=[(l['projection'],l['methods'][method]['observables'][observable]) for l in levels if method in l['methods']]
            if not rows:continue
            for value,v in rows:
                key={'mean':'mean','difference':'mean','sampling_difference':'mean','precision':'error','tau_int':'tau_int','ess':'ess','rhat':'rhat'}[quantity]
                y=v.get(key)
                if y is None:continue
                if quantity=='difference':
                    if v['ground_ed'] is None:continue
                    y-=v['ground_ed']
                if quantity=='sampling_difference':
                    if v['projected_ed'] is None:continue
                    y-=v['projected_ed']
                if quantity=='precision':y*=4
                error=v['error'] if quantity in ('mean','difference','sampling_difference') else v['tau_error'] if quantity=='tau_int' else None
                ax.errorbar(value,y,yerr=error,fmt='x' if v['flags'] else marker,color=color,capsize=3,markersize=5)
            ax.plot([],[],marker,color=color,label={'projector':'Open-boundary projector QMC','emc':'Independent classical EMC','sse':'Thermal SSE'}[method])
        if quantity=='precision' and levels[-1]['methods']:
            first=next(iter(levels[-1]['methods'].values()))['observables'][observable];ax.axhline(first['precision_tolerance'],color='black',linestyle='--',label='Required width')
        if quantity=='rhat':ax.axhline(1.01,color='black',linestyle='--',label='Review threshold')
        if quantity=='ess':ax.axhline(400,color='black',linestyle='--',label='Minimum ESS')
        if quantity in ('ess','tau_int'):ax.set_yscale('log')
    ax.set_xscale('log',base=2);ax.set_xlim(min(x)/1.2,max(x)*1.2)
    ax.set_xlabel(r'Total projection $T=2\tau$' if not point['kind'].startswith('classical') else r'Inverse temperature $\beta$')
    ax.set_ylabel({'mean':observable,'difference':'QMC - ground reference','sampling_difference':'QMC - exact finite ensemble','fidelity':'1 - exact projector ground-state fidelity','state_infidelity':'1 - reconstructed QMC-state fidelity','state_tvd':'Configuration total variation distance','state_missing':'Exact weight absent from QMC support','precision':'4 SE (review width)','tau_int':'Integrated autocorrelation (complete sweeps)','ess':'Effective sample size','rhat':'Split / folded R-hat'}[quantity])
    p=point['params'];ax.set_title(f'{p["Lx"]} x {p["Ly"]} periodic | '+', '.join(f'{k}={p[k]:.5g}' for k in ('J0','J1','J2','hx')),fontsize=9)
    ax.legend(frameon=False,fontsize=7)
    figure.suptitle(f'{point["label"]} | {point["verdict"]}',fontsize=9)
    figure.text(.02,.015,'QMC bars: 1 SE; crosses flag review. Finite-ensemble sampling and ground projection are separate checks.',fontsize=6)
    if quantity=='fidelity':figure.text(.02,.042,'Values below 1e-16 use the display floor; exact numerical values remain in JSON.',fontsize=6)
    if quantity.startswith('state_'):figure.text(.02,.042,'Distribution bars: independent-run bootstrap SE, conditional on observed support; not a fidelity guarantee.',fontsize=6)
    output=io.BytesIO();figure.savefig(output,format=format,dpi=300);return output.getvalue()
