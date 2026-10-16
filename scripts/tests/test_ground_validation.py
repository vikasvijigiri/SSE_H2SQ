# Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
from scipy.linalg import expm

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'backend'))
from app.services import ground_ed as ed
from app.services.ground_validation import Request,sample_parameters,evaluate_point,analytic_ground
from app.services.ground_statistics import summarize
from app.services.physics_model import basis_data
from app.services.physics_model import exact_reference


class GroundReferenceTests(unittest.TestCase):
    def test_projected_state_matches_independent_full_matrix_including_signed_couplings(self):
        for j0,j1,j2,h in [(1,.3,.025,.4),(-1,.3,.025,.4),(1,-.2,-.07,-.31)]:
            p=dict(Lx=2,Ly=2,J0=j0,J1=j1,J2=j2,hx=h)
            index,local,px,py,pairs,obs,_,_=basis_data(2,2)
            diagonal=-j0*np.prod(local,axis=2).sum(1)+j1*(local[:,:,0]*local[:,:,2]+local[:,:,1]*local[:,:,3]).sum(1)
            for a,b in pairs:diagonal+=j2*(px[:,a]*px[:,b]+py[:,a]*py[:,b])
            matrix=np.diag(diagonal);x=np.zeros_like(matrix)
            for i in range(4):x[index,index^(1<<i)]=1
            matrix-=h*x;e,v=np.linalg.eigh(matrix)
            trial=np.full(16,.25)
            if h<0:trial*=np.array([(-1)**int(i).bit_count() for i in index])
            for total in [.5,2.,8.]:
                wave=expm(-total/2*(matrix-e[0]*np.eye(16)))@trial;wave/=np.linalg.norm(wave)
                r=ed.reference(p,total)
                self.assertAlmostEqual(r['values']['E'],float(wave@matrix@wave)/4,places=12)
                self.assertAlmostEqual(r['values']['mx'],float(wave@x@wave)/4,places=12)
                for key in ['mz2','mz4','stripe2','stripe4','polar2','polar4','ice','plaquette','rho']:
                    self.assertAlmostEqual(r['values'][key],float(obs[key]@(wave*wave)),places=12,msg=key)
                self.assertGreaterEqual(r['values']['E'],e[0]/4-1e-12)
                self.assertAlmostEqual(r['ground_fidelity'],float(wave@v[:,0])**2,places=12)

    def test_midpoint_binomial_and_event_density_expansion_equal_exact_projector(self):
        index,_,_,_,_,obs,_,_=basis_data(2,2);p=dict(Lx=2,Ly=2,J0=-1.,J1=.3,J2=.025,hx=.4)
        diagonal=ed.diagonal_hamiltonian(2,2,-1.,.3,.025);h=np.diag(diagonal);x=np.zeros_like(h)
        for i in range(4):x[index,index^(1<<i)]=1
        h-=.4*x;shift=8.;k=shift*np.eye(16)-h;trial=np.full(16,.25);total=.5
        powers=[np.eye(16)]
        for i in range(40):powers.append(powers[-1]@k)
        z=0.;diag_sum=0.;event_sum=0.;count_sum=0.;a=np.diag(obs['mz2'])
        for n in range(41):
            factor=total**n/math.factorial(n);weight=factor*(trial@powers[n]@trial);z+=weight;count_sum+=n*weight
            for rank in range(n+1):diag_sum+=factor*math.comb(n,rank)/2**n*(trial@powers[n-rank]@a@powers[rank]@trial)
            for rank in range(n):event_sum+=factor*n/total*math.comb(n-1,rank)/2**(n-1)*(trial@powers[n-1-rank]@(.4*x)@powers[rank]@trial)
        wave=expm(total*k/2)@trial;norm=wave@wave
        self.assertAlmostEqual(z,norm,delta=norm*1e-11)
        self.assertAlmostEqual(diag_sum/z,float(wave@a@wave)/norm,places=11)
        self.assertAlmostEqual(event_sum/(.4*4*z),float(wave@x@wave)/(4*norm),places=11)
        self.assertAlmostEqual((shift-count_sum/(total*z))/4,float(wave@h@wave)/(4*norm),places=11)

    def test_random_points_frozen_balanced_and_simultaneously_varied(self):
        r=Request();a=sample_parameters(r);self.assertEqual(a,sample_parameters(r))
        self.assertEqual(len({tuple(np.sign(p[k]) for k in ['J1','J2','hx']) for p in a}),8)
        self.assertEqual(len({p['J1'] for p in a}),8);self.assertTrue(all(p['J1']*p['J2']*p['hx']!=0 for p in a))

    def test_wide_error_bars_never_earn_ground_go(self):
        r=Request(include_controls=False);reference=analytic_ground(4,4,.37)
        observations={key:dict(flags=[],precision_pass=False,projection_bias_pass=True,ground_agreement=True) for key in reference['values']}
        p=dict(kind='ground',levels=[dict(ed=reference,methods={'projector':dict(observables=observations)})])
        evaluate_point(p,r);self.assertEqual(p['verdict'],'NO GO');self.assertIn('requested statistical precision not reached',p['reasons'])

    def test_unavailable_exact_reference_remains_cross_check_only(self):
        rng=np.random.default_rng(27);traces=[]
        for i in range(4):
            a=np.empty(10000,dtype=[('E','f8')]);a['E']=rng.normal(size=len(a));traces.append(a)
        unknown={'values':{'E':None}}
        result=summarize(traces,unknown,unknown,1.)['observables']['E']
        self.assertIsNone(result['ground_ed']);self.assertIsNone(result['pull']);self.assertEqual(result['status'],'CROSS-CHECK ONLY')

    def test_one_converged_endpoint_is_insufficient_for_ground_go(self):
        r=Request(include_controls=False);reference=analytic_ground(4,4,.37)
        observations={key:dict(flags=[],precision_pass=True,projection_bias_pass=True,ground_agreement=True) for key in reference['values']}
        p=dict(kind='ground',levels=[dict(ed=reference,methods={'projector':dict(observables=observations)})])
        evaluate_point(p,r);self.assertEqual(p['verdict'],'NO GO');self.assertIn('projection-length plateau not measured',p['reasons'])

    def test_rotated_sixteen_spin_rectangles_have_identical_ground_observables(self):
        p=dict(J0=1.,J1=-.37,J2=.072,hx=.31)
        a=ed.reference(dict(Lx=2,Ly=8,**p));b=ed.reference(dict(Lx=8,Ly=2,**p))
        for key in ed.KEYS:self.assertAlmostEqual(a['values'][key],b['values'][key],delta=3e-11,msg=key)
        self.assertLess(a['numerical_residual'],1e-11);self.assertLess(b['numerical_residual'],1e-11)

    @unittest.skipUnless(os.environ.get('SSE_FULL_REFERENCE_TESTS')=='1','Full 16-spin ground check')
    def test_full_4x4_ground_matches_analytic_product_limit(self):
        p=dict(Lx=4,Ly=4,J0=0.,J1=0.,J2=0.,hx=-.37)
        expected=analytic_ground(4,4,-.37);actual=ed.reference(p,2.)
        self.assertEqual(actual['dimension'],65536)
        for key,value in expected['values'].items():self.assertAlmostEqual(value,actual['values'][key],places=11,msg=key)


class GroundEngineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(['make','-C',str(ROOT/'backend/build'),'projector','benchmark'],check=True,capture_output=True)

    def test_conditional_estimators_equal_explicit_component_and_symmetry_averages(self):
        source=r'''
#define main projector_runner_main
#include "backend/benchmarks/projector.cpp"
#undef main
int main(){
 std::mt19937 random(61947);
 for(auto shape:std::vector<std::array<int,2>>{{2,2},{4,2}}){
  Parameters p;p.Lx=shape[0];p.Ly=shape[1];p.Lz=1;p.Ns=p.Lx*p.Ly;p.J0=-1;p.J1=.31;p.J2=-.071;
  Geometry g;g.build(p.Lx,p.Ly,1);ProjectorDiagonal estimator(p,g);
  std::vector<int> mapping(1<<p.Ns);std::iota(mapping.begin(),mapping.end(),0);
  std::vector<double> histogram(mapping.size());ProjectorHistogram hist(mapping,histogram);
  for(int sample=0;sample<200;++sample){
   int k=1+random()%p.Ns;std::vector<int> spin(p.Ns),component(p.Ns);for(int i=0;i<p.Ns;++i){spin[i]=random()%2?1:-1;component[i]=i<k?i:random()%k;}
   std::array<double,11> expected{};std::vector<double> expected_hist(mapping.size());
   for(int mask=0;mask<(1<<k);++mask){auto s=spin;for(int i=0;i<p.Ns;++i)if(mask&(1<<component[i]))s[i]*=-1;unsigned state=0;for(int i=0;i<p.Ns;++i)if(s[i]>0)state|=1u<<i;expected_hist[state]+=1./(1<<k);
    for(int orbit=0;orbit<4;++orbit){auto t=s;for(int i=0;i<p.Ns;++i){int parity=((i%p.Lx)+(i/p.Lx))%2;if(orbit&(1<<parity))t[i]*=-1;}auto v=diagonal_measure(t,p,g);for(int j=0;j<11;++j)expected[j]+=v[j]/(4*(1<<k));}
   }
   auto actual=estimator.measure(spin,component);for(int j=0;j<11;++j)if(std::abs(expected[j]-actual[j])>2e-12){std::cerr<<"conditional mismatch "<<p.Ns<<" "<<sample<<" "<<j<<" "<<expected[j]<<" "<<actual[j];return 1;}
   std::fill(histogram.begin(),histogram.end(),0.);hist.add(spin,component);for(size_t j=0;j<histogram.size();++j)if(std::abs(histogram[j]-expected_hist[j])>1e-14)return 2;
  }
 }
 return 0;
}
'''
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp);(folder/'contract.cpp').write_text(source)
            subprocess.run(['g++','-O2','-std=c++17','-I',str(ROOT),str(folder/'contract.cpp'),str(ROOT/'backend/src/updates.cpp'),str(ROOT/'backend/src/variables.cpp'),str(ROOT/'backend/src/observables.cpp'),'-static-libgcc','-static-libstdc++','-o',str(folder/'contract.exe')],check=True,capture_output=True)
            result=subprocess.run([str(folder/'contract.exe')],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)

    def test_independent_classical_energy_for_every_basis_state(self):
        with tempfile.TemporaryDirectory() as tmp:
            f=Path(tmp)
            for lx,ly,j0,j1,j2 in [(2,2,-1.,.3,.025),(4,2,1.,-.23,.071),(4,4,1.,.37,-.043)]:
                p=dict(Lx=lx,Ly=ly,J0=j0,J1=j1,J2=j2,hx=0.,beta=1.,Lz=1,J3=0,n_therm=0,n_measure=1)
                (f/'input').write_text('\n'.join(f'{k} {v}' for k,v in p.items()))
                subprocess.run([str(ROOT/'backend/build/emc_validation.exe'),str(f/'input'),str(f/'trace'),'34','enumerate'],check=True,capture_output=True)
                a=pd.read_csv(f/'trace')
                np.testing.assert_allclose(a['E']*(lx*ly),ed.diagonal_hamiltonian(lx,ly,j0,j1,j2),atol=2e-13,rtol=0)

    def test_open_projector_with_frustration_matches_finite_projected_ed(self):
        p=dict(Lx=2,Ly=2,J0=-1.,J1=.3,J2=.025,hx=.4)
        with tempfile.TemporaryDirectory() as tmp:
            f=Path(tmp);(f/'input').write_text('\n'.join(f'{k} {v}' for k,v in {**p,'beta':2.,'Lz':1,'J3':0,'n_therm':10000,'n_measure':100000}.items()))
            mapping,probabilities=ed.orbit_probabilities(p,2.);mapping.astype('<i4').tofile(f/'map')
            def run(i):
                import json
                result=subprocess.run([str(ROOT/'backend/build/sse_projector.exe'),str(f/'input'),str(f/f'trace-{i}'),str(1743+i*1009),str(f/'map')],check=True,capture_output=True,text=True)
                metadata=json.loads(result.stdout);self.assertGreater(metadata['open_endpoint_rows'],1000);self.assertEqual(metadata['cutoff_saturations'],0)
                counts=np.fromfile(f/f'trace-{i}.hist.bin',dtype='<f8');self.assertAlmostEqual(counts.sum(),100000,places=6)
                self.assertGreater(np.sum(np.sqrt(counts/counts.sum()*probabilities))**2,.999)
                frame=pd.read_csv(f/f'trace-{i}');return frame.to_records(index=False)
            with ThreadPoolExecutor(max_workers=4) as pool:traces=list(pool.map(run,range(4)))
            summary=summarize(traces,ed.reference(p,2.),ed.reference(p),30.)
            for key,v in summary['observables'].items():self.assertFalse(v['flags'],(key,v))

    def test_independent_exchange_mc_and_thermal_sse_both_match_exact_classical_trace(self):
        import json
        p=dict(Lx=4,Ly=2,J0=1.,J1=-.23,J2=.071,hx=0.,beta=.5)
        ref=exact_reference(p);ref['values']['E_count']=ref['values']['E']
        ref['values'].update({key:ref['values'][original] for key,original in ed.RAW_ALIASES.items()})
        with tempfile.TemporaryDirectory() as tmp:
            f=Path(tmp);(f/'input').write_text('\n'.join(f'{k} {v}' for k,v in {**p,'Lz':1,'J3':0,'n_therm':2000,'n_measure':10000}.items()))
            for engine in ['emc_validation','sse_benchmark']:
                def run(i):
                    command=[str(ROOT/f'backend/build/{engine}.exe'),str(f/'input'),str(f/f'{engine}-{i}')]
                    command += [str(1523+i*1009),'run'] if engine=='emc_validation' else ['symmetry',str(2519+i*1009)]
                    r=subprocess.run(command,check=True,capture_output=True,text=True);meta=json.loads(r.stdout)
                    if engine=='emc_validation':self.assertGreater(meta['round_trips'],10)
                    frame=pd.read_csv(f/f'{engine}-{i}');a=np.empty(len(frame),dtype=[(key,'f8') for key in ed.KEYS])
                    for key in ed.KEYS:a[key]=frame[key if key in frame else 'E' if key=='E_count' else ed.RAW_ALIASES[key]]
                    return a
                with ThreadPoolExecutor(max_workers=4) as pool:traces=list(pool.map(run,range(4)))
                result=summarize(traces,ref,ref,1.,analytic_constants=['mx','raw_mx'])
                for key,v in result['observables'].items():self.assertFalse(v['flags'],(engine,key,v))


if __name__=='__main__':unittest.main()
