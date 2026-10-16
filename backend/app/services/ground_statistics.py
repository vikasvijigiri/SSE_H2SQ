# Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
"""Independent-run precision and pure-state convergence criteria."""
import math
import numpy as np
from .benchmark_statistics import autocorrelation,split_rhat


def summarize(traces,reference,ground,total_seconds,precision=0.01,analytic_constants=()):
    sizes=np.asarray([len(t) for t in traces],dtype=float);weights=sizes/sizes.sum();n=int(sizes.sum())
    observations={};streams={};largest=.5
    for key in traces[0].dtype.names:
        stream=[t[key] for t in traces];stats=[autocorrelation(s) for s in stream]
        taus=[s['tau'] for s in stats if s['tau'] is not None]
        largest=max(largest,max(taus,default=.5));streams[key]=(stream,stats)
    requested=max(32,math.ceil(10*largest));block=min(requested,min(map(len,traces))//8)
    for key,(stream,stats) in streams.items():
        means=np.asarray([s.mean() for s in stream]);mean=float(weights@means)
        blocks=np.concatenate([s[:len(s)//block*block].reshape(-1,block).mean(1) for s in stream])
        block_se=float(np.std(blocks,ddof=1)*np.sqrt(block/n))
        between=math.sqrt(float(np.sum(weights**2*(means-mean)**2/(1-weights))))
        acf_se=math.sqrt(sum(len(s)*np.var(s,ddof=1)*2*r['tau'] for s,r in zip(stream,stats) if r['tau'] is not None)/n**2)
        se=max(block_se,between,acf_se);rhat=split_rhat(stream);ess=sum(r['ess'] or 0 for r in stats)
        constant=all(r['constant'] for r in stats);flags=[]
        target=reference['values'][key];gs=ground['values'][key]
        # Target absolute precision complements residual checks: broad bars do not earn GO.
        floors={'E':1e-3,'E_count':1e-3,'E2':2e-3,'E4':4e-3,'mz2':1e-5,'mz4':1e-7}
        tolerance=max(floors.get(key,1e-3),precision*abs(gs or 0.))
        if rhat is not None and (not math.isfinite(rhat) or rhat>1.01):flags.append('between-run disagreement')
        analytic=constant and key in analytic_constants
        if constant and not analytic:flags.append('constant trace: mixing unresolved')
        if ess<400 and not constant:flags.append('low effective sample size')
        if requested>min(map(len,traces))//8:flags.append('insufficient correlation-scale blocks')
        if len(traces)<4:flags.append('fewer than four independent runs')
        residual=mean-target if target is not None else None;pull=residual/se if residual is not None and se>0 else None
        if residual is not None and abs(residual)>4*se+1e-10:flags.append('finite-ensemble exact-reference discrepancy exceeds 4 SE')
        bias=abs(target-gs) if target is not None and gs is not None else None
        ground_agreement=abs(mean-gs)<=4*se+1e-10 if gs is not None else None
        precision_pass=(not constant or analytic) and 4*se<=tolerance
        bias_pass=bias<=tolerance/4 if bias is not None else None
        taus=[r['tau'] for r in stats if r['tau'] is not None]
        observations[key]=dict(mean=mean,error=se if not constant or analytic else None,tau_int=float(np.mean(taus)) if taus else None,
            tau_error=float(np.mean([r['tau_error'] for r in stats if r['tau_error'] is not None])) if taus else None,
            ess=ess,rhat=rhat if rhat is not None and math.isfinite(rhat) else None,flags=flags,pull=pull,
            projected_ed=target,ground_ed=gs,projection_bias=bias,precision_tolerance=tolerance,
            precision_pass=precision_pass,projection_bias_pass=bias_pass,ground_agreement=ground_agreement,
            status='NO GO' if flags else 'CROSS-CHECK ONLY' if target is None else 'GO' if precision_pass and bias_pass and ground_agreement else 'NOT CONVERGED',
            ess_per_second=ess/total_seconds if total_seconds else None,
            acf=[dict(lag=i,value=float(np.mean([r['acf'][i] for r in stats if r['acf'] is not None]))) for i in [0,1,2,4,8,16,32,64,128,256] if i<min(map(len,traces))] if taus else [])
    return dict(observables=observations,retained_samples=n,chains=len(traces),block_size=block,
        requested_block_size=requested,total_seconds=total_seconds,mean_method='sample-count weighted; all rows retained',
        error_method='maximum of blocked, chain-wise autocorrelation and independent-run HC2 SE')
