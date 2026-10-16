# Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
"""Independent checks of the complete-sector reference and improved estimators."""
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"backend"))
from app.services.physics_model import basis_data,exact_reference
from app.services.symmetry_ed import sectors,hamiltonian,exact_4x4_reference
from app.services.benchmark_statistics import PRIMARY,summarize


class SymmetryReferenceTests(unittest.TestCase):
    def test_all_sector_dimensions_and_spectral_trace_moments(self):
        blocks,_,group=sectors()
        indices,local,px,py,pairs,_,_,_=basis_data(4,4)
        diagonal=-np.prod(local,axis=2).sum(1)-.23*(local[:,:,0]*local[:,:,2]+local[:,:,1]*local[:,:,3]).sum(1)
        for a,b in pairs:diagonal+=.047*(px[:,a]*px[:,b]+py[:,a]*py[:,b])
        h,_=hamiltonian(diagonal,.37,16)
        self.assertEqual(sum(s['basis'].shape[1] for s in blocks),65536)
        trace,second=0.,0.
        for item in blocks:
            b=item['basis'];projected=b.T@h@b
            np.testing.assert_allclose((b.T@b).diagonal(),1.,atol=1e-14)
            trace+=projected.diagonal().sum();second+=projected.multiply(projected).sum()
        self.assertAlmostEqual(trace,diagonal.sum(),places=8)
        self.assertAlmostEqual(second,float(diagonal@diagonal)+65536*16*.37**2,places=7)
        for p in group[[1,2,4,8,16]]:np.testing.assert_allclose(diagonal,diagonal[p],atol=3e-14)

    def test_lehmann_responses_and_static_moments_match_independent_dense_ed(self):
        for hx,j1,j2,beta in [(.2,.3,.025,2.),(-.31,-.2,.04,3.),(.4,.2,-.03,80.)]:
            p=dict(Lx=2,Ly=2,J0=1.,J1=j1,J2=j2,hx=hx,beta=beta)
            dense=exact_reference(p);sector=exact_4x4_reference(p)
            self.assertEqual(sector['dimension'],16)
            for key,value in dense['values'].items():
                if value is not None:
                    bound=sector['absolute_error_bounds'].get(key,0.)
                    self.assertAlmostEqual(value,sector['values'][key],delta=bound+3e-11,msg=key)

    @unittest.skipUnless(os.environ.get('SSE_FULL_REFERENCE_TESTS')=='1','Set SSE_FULL_REFERENCE_TESTS=1 for the full 16-spin analytic check')
    def test_full_4x4_noninteracting_transverse_field(self):
        h,beta=.37,32.
        p=dict(Lx=4,Ly=4,J0=0.,J1=0.,J2=0.,hx=h,beta=beta)
        result=exact_reference(p);v=result['values'];mx=np.tanh(beta*h)
        self.assertEqual(result['dimension'],65536)
        self.assertAlmostEqual(v['E'],-h*mx,places=11)
        self.assertAlmostEqual(v['mx'],mx,places=11)
        self.assertAlmostEqual(v['mz2'],1/64,places=11)
        self.assertAlmostEqual(v['stripe2'],1/32,places=11)
        self.assertAlmostEqual(v['chi_stag'],mx/(4*h),places=10)
        self.assertAlmostEqual(v['chi_x'],beta/np.cosh(beta*h)**2,delta=2e-10)
        for i,f in enumerate([0,.125,.25,.375,.5]):
            expected=.25*np.cosh(beta*h*(1-2*f))/np.cosh(beta*h)
            self.assertAlmostEqual(v[f'ct{i}'],expected,delta=result['absolute_error_bounds'].get(f'ct{i}',0.)+3e-11)


class PrecisionTests(unittest.TestCase):
    def traces(self,lengths):
        columns=list(dict.fromkeys(PRIMARY+['stripe4','bar_m','bar_m2','bar_px','bar_py','bar_p2','nx','nxf2']))
        rng=np.random.default_rng(9381);traces=[]
        for n in lengths:
            trace=np.zeros(n,dtype=[(key,'f8') for key in columns])
            trace['E']=rng.normal(size=n);trace['E2']=trace['E']**2;trace['E4']=trace['E']**4
            trace['stripe2']=.1;trace['stripe4']=.01
            traces.append(trace)
        return traces

    def test_unequal_chain_lengths_use_every_row(self):
        traces=self.traces([2000,4000,3000,5000])
        traces[1]['E']+=.5
        p=dict(Lx=2,Ly=2,beta=1.,hx=0.)
        result=summarize(traces,p)
        expected=sum(t['E'].sum() for t in traces)/14000
        self.assertAlmostEqual(result['observables']['E']['mean'],expected,places=14)
        self.assertEqual(result['retained_samples'],14000)
        self.assertEqual(result['discarded_tail_per_chain'],0)

    def test_independent_gaussian_se_scales_with_sample_count(self):
        traces=self.traces([12000]*12);p=dict(Lx=2,Ly=2,beta=1.,hx=0.)
        first=summarize(traces[:4],p)['observables']['E']['error']
        last=summarize(traces,p)['observables']['E']['error']
        self.assertLess(last,first)
        # The max-of-estimates rule can inflate a four-chain between-run SE;
        # it must not be forced to follow an exact 1/sqrt(N) trajectory.
        self.assertAlmostEqual(last,1/np.sqrt(144000),delta=.25/np.sqrt(144000))

    @unittest.skipUnless(shutil.which('g++'),'C++ compiler unavailable')
    def test_piecewise_path_estimator_matches_explicit_cyclic_binomial_average(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'.cache') as folder:
            exe=Path(folder)/('estimator.exe' if os.name=='nt' else 'estimator')
            command=['g++','-O2','-std=c++17','backend/benchmarks/estimator_contract.cpp','-o',str(exe)]
            if os.name=='nt':command+=['-static-libgcc','-static-libstdc++']
            subprocess.run(command,cwd=ROOT,check=True,capture_output=True)
            error=float(subprocess.check_output([str(exe)],text=True))
            self.assertLess(error,2e-12)


class PoolingTests(unittest.TestCase):
    def fixture(self,identity,seeds,lengths,method='canonical'):
        from app.services import benchmark
        params=dict(Lx=2,Ly=2,J0=1.,J1=.3,J2=.025,hx=0.,beta=2.)
        traces=PrecisionTests().traces(lengths)
        runs=[]
        for chain,(seed,trace) in enumerate(zip(seeds,traces)):
            folder=benchmark.STORE/identity/('a'*16)/method/str(chain)
            folder.mkdir(parents=True);np.savez_compressed(folder/'trace.npz',trace=trace)
            runs.append(dict(seed=seed,wall_seconds=1.,measurement_seconds=1.,cutoff_saturations=0,update_seconds=dict(diagonal=0.,cluster=0.,gauge=0.,canonical=0.)))
        job=dict(id=identity,created='2026-10-16',status='complete',request={**params,'axes':{'beta':[2.]},'algorithms':[method],'chains':len(seeds),'n_measure':lengths[0],'n_therm':100,'seed':seeds[0]},completed_runs=len(seeds),total_runs=len(seeds),points=[dict(axis='beta',x=2.,params=params,point_id='a'*16,ed={'values':{'E':0.}},methods={method:{'runs':runs}})],provenance={'model_version':'test','source_sha256':{}},definitions={},citations=[])
        benchmark.save_job(job)
        return traces

    def test_pool_recomputes_all_rows_and_rejects_repeated_seeds_and_aliases(self):
        from app.services import benchmark
        with tempfile.TemporaryDirectory(dir=ROOT/'.cache') as folder,patch.object(benchmark,'STORE',Path(folder)),patch.object(benchmark,'provenance',return_value={'model_version':'test','source_sha256':{}}),patch.object(benchmark,'archive_sources'):
            a=self.fixture('1'*32,[11,12],[1200,1600]);b=self.fixture('2'*32,[21,22],[1800,2400])
            pooled=benchmark.pool_jobs(['1'*32,'2'*32]);summary=pooled['points'][0]['methods']['canonical']
            self.assertEqual(summary['retained_samples'],7000)
            self.assertAlmostEqual(summary['observables']['E']['mean'],sum(t['E'].sum() for t in a+b)/7000,places=14)
            self.assertEqual(benchmark.trace_location(pooled['id'],pooled['points'][0],'canonical',3),('2'*32,'a'*16,1))
            source=benchmark.get_job('2'*32);source['points'][0]['methods']['canonical']['runs'][0]['seed']=11;benchmark.save_job(source)
            with self.assertRaisesRegex(RuntimeError,'Repeated seed'):benchmark.pool_jobs(['1'*32,'2'*32])
            source['points'][0]['methods']['canonical']['runs'][0].update(source_experiment='1'*32,source_point='a'*16,source_chain=0);benchmark.save_job(source)
            with self.assertRaisesRegex(RuntimeError,'selected twice'):benchmark.pool_jobs(['1'*32,'2'*32])

    def test_completed_primary_retains_stopped_controls_without_fabricating_completion(self):
        from app.services import benchmark
        with tempfile.TemporaryDirectory(dir=ROOT/'.cache') as folder,patch.object(benchmark,'STORE',Path(folder)):
            self.fixture('1'*32,[11,12],[1200,1200],method='symmetry')
            self.fixture('2'*32,[21,22],[1200,1200])
            control=benchmark.get_job('2'*32);control.update(status='cancelled',total_runs=8);benchmark.save_job(control)
            combined=benchmark.overlay_controls('1'*32,['2'*32])
            self.assertEqual(combined['completed_runs'],4)
            self.assertEqual(combined['source_experiments'][1]['status'],'cancelled')
            self.assertEqual(benchmark.get_job('2'*32)['status'],'cancelled')
            self.assertEqual(set(combined['points'][0]['methods']),{'canonical','symmetry'})
            self.assertEqual(benchmark.trace_location(combined['id'],combined['points'][0],'canonical',1),('2'*32,'a'*16,1))
            control['points'][0]['ed']['values']['E']=1.;benchmark.save_job(control)
            with self.assertRaisesRegex(RuntimeError,'ED references disagree'):benchmark.overlay_controls('1'*32,['2'*32])


@unittest.skipUnless(shutil.which('g++'),'C++ compiler unavailable')
class SiteSymmetryTests(unittest.TestCase):
    def test_every_shifted_diagonal_weight_and_symmetry_breaking_guard(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'.cache') as folder:
            exe=Path(folder)/('symmetry.exe' if os.name=='nt' else 'symmetry')
            command=['g++','-O2','-std=c++17','backend/benchmarks/symmetry_contract.cpp','backend/src/updates.cpp','backend/src/variables.cpp','backend/src/observables.cpp','-o',str(exe)]
            if os.name=='nt':command+=['-static-libgcc','-static-libstdc++','-Wl,-Bstatic,--whole-archive','-lwinpthread','-Wl,--no-whole-archive,-Bdynamic']
            subprocess.run(command,cwd=ROOT,check=True,capture_output=True)
            for lx,ly in [(2,2),(4,2),(2,4),(4,4)]:
                self.assertEqual(float(subprocess.check_output([str(exe),str(lx),str(ly)],text=True)),0.)

    def test_quantum_symmetry_schedule_closes_worldlines_and_matches_dense_ed(self):
        from app.services.benchmark import BenchmarkRequest,run_chain,EXECUTABLE
        import threading
        if not EXECUTABLE.exists():self.skipTest('Build benchmark engine first')
        request=BenchmarkRequest(Lx=4,Ly=2,beta=2.,J1=-.3,J2=-.025,hx=.3,algorithms=['symmetry'],chains=4,n_therm=2000,n_measure=12000,axes={'beta':[2.]})
        params={k:getattr(request,k) for k in ['Lx','Ly','J0','J1','J2','hx','beta']}
        with tempfile.TemporaryDirectory(dir=ROOT/'.cache') as folder:
            runs=[run_chain(params,'symmetry',77100+i,Path(folder)/str(i),request,threading.Event()) for i in range(4)]
            for _,metadata in runs:
                self.assertEqual(sum(metadata['site_symmetry_orbit_counts']),12000)
                self.assertTrue(all(count>0 for count in metadata['site_symmetry_orbit_counts']))
                self.assertEqual(metadata['cutoff_saturations'],0)
                self.assertGreaterEqual(metadata['configuration_checks'],14)
            result=summarize([trace for trace,_ in runs],params)
        reference=exact_reference(params)
        for key in ['E','mz2','ice','stripe2','polar2','ct0','ct4','chi_stag','chi_polar']:
            observed=result['observables'][key]
            self.assertLess(abs(observed['mean']-reference['values'][key]),8*observed['error']+2e-4,msg=key)


if __name__=='__main__':unittest.main()
