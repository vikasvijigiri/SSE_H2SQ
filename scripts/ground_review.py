# Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
"""Independently audit retained ground-study rows, hashes and publication assets."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('identity');parser.add_argument('--paper',type=Path)
    args=parser.parse_args();folder=ROOT/'.cache/ground-validation'/args.identity
    job=json.loads((folder/'study.json').read_text());assert job['status']=='complete'
    for name,sha in job['provenance']['source_sha256'].items():assert digest(folder/'sources'/sha/name)==sha,name
    for name,sha in job['provenance']['engine_sha256'].items():assert digest(folder/'engines'/f'{name}.exe')==sha,name
    seeds=set();runs=0;rows=0;pulls=[];rhats=[];taus=[];effective=[];ratios=[];flags=[];distributions=[]
    for point in job['points']:
        if point.get('orbit_reference'):
            reference=point['orbit_reference']
            assert digest(folder/reference['map_file'])==reference['map_sha256']
            assert digest(folder/reference['ground_file'])==reference['ground_sha256']
        for level in point['levels']:
            if level.get('orbit_probability_file'):assert digest(folder/level['orbit_probability_file'])==level['orbit_probability_sha256']
            for method,summary in level['methods'].items():
                totals={key:0. for key in summary['observables']};n=0
                for index,metadata in enumerate(summary['runs']):
                    assert metadata['seed'] not in seeds,'Repeated seed';seeds.add(metadata['seed']);runs+=1
                    path=folder/point['id']/f'{level["projection"]:g}'/method/str(index)
                    stored=json.loads((path/'metadata.json').read_text());assert all(metadata.get(key)==value for key,value in stored.items())
                    inputs=dict(line.split(maxsplit=1) for line in (path/'input.in').read_text().splitlines())
                    assert int(inputs['n_measure'])==metadata['measurement_rows'];assert float(inputs['beta'])==level['projection']
                    for key,value in point['params'].items():
                        if key!='beta':assert float(inputs[key])==value,(point['id'],key)
                    with np.load(path/'trace.npz',allow_pickle=False) as saved:
                        trace=saved['trace'];assert len(trace)==metadata['measurement_rows'];rows+=len(trace);n+=len(trace)
                        for key in totals:assert np.isfinite(trace[key]).all();totals[key]+=float(np.sum(trace[key]))
                        if method=='projector' and point.get('orbit_reference'):
                            histogram=saved['histogram'];assert np.all(histogram>=0);assert abs(histogram.sum()-len(trace))<1e-6
                    if method!='emc':assert metadata['cutoff_saturations']==0
                assert n==summary['retained_samples']
                for key,v in summary['observables'].items():
                    assert abs(totals[key]/n-v['mean'])<1e-11*max(1.,abs(v['mean'])),(point['id'],method,key)
                    if v['pull'] is not None:pulls.append(abs(v['pull']))
                    if v['rhat'] is not None:rhats.append(v['rhat'])
                    if v['tau_int'] is not None:taus.append(v['tau_int']);effective.append(v['ess'])
                    if v['flags']:flags.append(dict(point=point['id'],T=level['projection'],method=method,observable=key,flags=v['flags']))
                if summary.get('precision_history'):
                    ratios.extend(v['SE_ratio'] for v in summary['precision_history'][-1]['observables'].values() if v['SE_ratio'] is not None)
                if summary.get('histogram'):
                    h=summary['histogram'];assert abs(h['fidelity_ground']-h['reconstructed_state']['direct_ground_overlap'])<1e-10
                    distributions.append(dict(point=point['id'],T=level['projection'],**{k:h[k] for k in ['fidelity_ground','fidelity_projected','tvd_ground','missing_ground_weight']}))
    assert runs==job['completed_runs']==job['total_runs'];assert len(seeds)==runs
    assert job['verdict']['universal']=='NO GO'
    requested=[p for p in job['points'] if p['label'].startswith(('Random','Custom'))]
    assert job['verdict']['random_ground_domain']==('GO' if requested and all(p['verdict']=='GO' for p in requested) else 'NO GO')
    report=dict(id=job['id'],source_revision=job['provenance']['revision'],runs=runs,rows=rows,unique_seeds=len(seeds),
        universal_verdict=job['verdict']['universal'],requested_domain=job['verdict']['random_ground_domain'],
        maximum_absolute_ED_pull=max(pulls,default=None),maximum_rhat=max(rhats,default=None),minimum_ess=min(effective,default=None),
        tau_range=[min(taus,default=None),max(taus,default=None)],SE_ratios=dict(count=len(ratios),decreased=sum(r<1 for r in ratios),median=float(np.median(ratios)) if ratios else None),
        cases=[dict(id=p['id'],label=p['label'],verdict=p['verdict'],reasons=p.get('reasons',[])) for p in job['points']],flags=flags,distributions=distributions)
    if args.paper:
        paper=args.paper.resolve();assert paper.is_relative_to(ROOT)
        manifest=json.loads((paper/'metadata/manifest.json').read_text())
        for name,sha in manifest.items():assert digest(paper/name)==sha,name
        for source in paper.rglob('*.tex'):
            text=re.sub(r'(?<!\\)%[^\n]*','',source.read_text(encoding='utf-8'))
            for pattern,extension in [(r'\\input\{([^}]+)\}','.tex'),(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}',''),(r'\\bibliography\{([^}]+)\}','.bib')]:
                for name in re.findall(pattern,text):assert (paper/(name if Path(name).suffix else name+extension)).is_file(),name
        bib='\n'.join(p.read_text() for p in paper.rglob('*.bib'));keys=set(re.findall(r'@\w+\s*\{\s*([^,\s]+)',bib))
        citations=set(k.strip() for source in paper.rglob('*.tex') for group in re.findall(r'\\cite\w*\{([^}]+)\}',source.read_text()) for k in group.split(','));assert citations<=keys
        for path in (paper/'figures').glob('*'):
            data=path.read_bytes()
            assert (path.suffix=='.pdf' and data.startswith(b'%PDF-')) or (path.suffix=='.png' and data.startswith(b'\x89PNG')) or (path.suffix=='.svg' and b'<svg' in data[:1000]),path
        report['publication_assets_checked']=len(manifest);report['latex_compiled']=False
        if (paper/'metadata/latex-build.json').exists():
            build=json.loads((paper/'metadata/latex-build.json').read_text());assert build['compiled'];assert digest(paper/'validation-report.pdf')==build['pdf_sha256']
            assert build.get('input_sha256'),'Compiled inputs are not recorded.'
            assert all(digest(paper/name)==sha for name,sha in build['input_sha256'].items()),'Compiled PDF is stale.'
            report['latex_compiled']=True;report['latex_pages']=build['pages']
        (paper/'metadata/review.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
        retained={p.relative_to(folder).as_posix():digest(p) for p in sorted(folder.rglob('*')) if p.is_file() and p.name!='study.tmp'}
        (paper/'metadata/retained-artifacts.json').write_text(json.dumps(retained,indent=2),encoding='utf-8')
        manifest={p.relative_to(paper).as_posix():digest(p) for p in sorted(paper.rglob('*')) if p.is_file() and p.name!='manifest.json'}
        (paper/'metadata/manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    target=ROOT/'.cache'/f'ground-review-{args.identity}.json';target.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('flags','distributions','cases')},indent=2))


if __name__=='__main__':main()
