# Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
"""Pure-state ground and finite-projector references; never a thermal trace."""
from functools import lru_cache
import hashlib
import inspect
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from .physics_model import basis_data, geometry
from .symmetry_ed import hamiltonian, _permutation

ROOT=Path(__file__).resolve().parents[3]
CACHE=ROOT/'.cache/ground-ed'
KEYS=['E','E_count','E2','E4','mz2','mz4','mx','plaquette','ice','rho','stripe2','stripe4','polar2','polar4','ct0']
RAW_ALIASES={'raw_mx':'mx','raw_mz2':'mz2','raw_mz4':'mz4','raw_ice':'ice'}
KEYS+=list(RAW_ALIASES)


@lru_cache(maxsize=8)
def positive_basis(lx,ly):
    """Invariant positive orbit basis. PF fixes the unique h!=0 ground sector."""
    n=lx*ly
    if lx%2 or ly%2 or min(lx,ly)<2 or n>16:
        raise ValueError('Ground ED supports even periodic cells with at most 16 spins.')
    indices=np.arange(1<<n,dtype=np.int32)
    even=sum(1<<i for i in range(n) if (i%lx+i//lx)%2==0)
    # Ground/trial invariance needs no commuting-character construction.
    # The full orbit group works for rectangular cells, including translations
    # of order three on 2x6. Verify its action on the actual Hamiltonian below.
    spatial=[]
    for dx in range(0,lx,2):
        for dy in range(0,ly,2):
            for reflect_x in (False,True):
                for reflect_y in (False,True):
                    mapping=[]
                    for site in range(n):
                        x,y=site%lx,site//lx
                        tx=((1-x if reflect_x else x)+dx)%lx
                        ty=((1-y if reflect_y else y)+dy)%ly
                        mapping.append(tx+lx*ty)
                    spatial.append(_permutation(indices,mapping))
    spatial=np.asarray(spatial)
    if lx==ly:
        transpose=_permutation(indices,[site//lx+lx*(site%lx) for site in range(n)])
        spatial=np.concatenate([spatial,spatial[:,transpose]])
    group=np.concatenate([spatial^mask for mask in (0,even,((1<<n)-1)^even,(1<<n)-1)])
    reps=np.min(group,axis=0);unique,inverse,counts=np.unique(reps,return_inverse=True,return_counts=True)
    basis=sparse.csc_matrix((1/np.sqrt(counts[inverse]),(indices,inverse)),shape=(len(indices),len(unique)))
    return basis,group


def diagonal_hamiltonian(lx,ly,j0,j1,j2):
    _,local,px,py,pairs,_,_,_=basis_data(lx,ly)
    diagonal=-j0*np.prod(local,axis=2).sum(1)+j1*(local[:,:,0]*local[:,:,2]+local[:,:,1]*local[:,:,3]).sum(1)
    for a,b in pairs:diagonal+=j2*(px[:,a]*px[:,b]+py[:,a]*py[:,b])
    return diagonal


@lru_cache(maxsize=2)
def prepared(lx,ly,j0,j1,j2,habs):
    if habs<=0:raise ValueError('Unique positive ground-state validation requires nonzero transverse field.')
    basis,group=positive_basis(lx,ly)
    if basis.shape[1]>4096:
        raise ValueError('Exact projector block exceeds the dense memory budget; this point is not ED validated.')
    diagonal=diagonal_hamiltonian(lx,ly,j0,j1,j2)
    if any(not np.allclose(diagonal,diagonal[g],atol=1e-12,rtol=0) for g in group):
        raise RuntimeError('Ground-sector projection does not preserve the actual Hamiltonian.')
    h,transverse=hamiltonian(diagonal,habs,lx*ly)
    block=(basis.T@h@basis).toarray()
    fingerprint=hashlib.sha256((str((lx,ly,j0,j1,j2,habs))+inspect.getsource(positive_basis.__wrapped__)+inspect.getsource(diagonal_hamiltonian)).encode()+block.tobytes()).hexdigest()[:24]
    path=CACHE/f'{fingerprint}.npz';path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():
        with np.load(path,allow_pickle=False) as z:energies=z['energies'];vectors=z['vectors']
    else:
        energies,vectors=eigh(block,driver='evd',check_finite=False)
        np.savez(path,energies=energies,vectors=vectors)
    residual=float(np.max(np.linalg.norm(block@vectors-vectors*energies,axis=0)))
    ground=basis@vectors[:,0]
    residual_full=float(np.linalg.norm(h@ground-energies[0]*ground))
    if max(residual,residual_full)>1e-9:raise RuntimeError('Ground ED residual is too large.')
    if np.sum(ground)<0:ground=-ground
    if ground.min() < -1e-10:raise RuntimeError('Lowest invariant state fails the positive-ground-state check.')
    trial=np.asarray(basis.T@np.full(1<<(lx*ly),1/np.sqrt(1<<(lx*ly))))
    coefficients=vectors.T@trial
    return basis,h,transverse,diagonal,energies,vectors,coefficients,max(residual,residual_full)


def reference(params,total_projection=None):
    lx,ly=int(params['Lx']),int(params['Ly']);n=lx*ly;hx=float(params['hx'])
    b,h,transverse,diagonal,e,v,c,residual=prepared(lx,ly,float(params['J0']),float(params['J1']),float(params['J2']),abs(hx))
    if total_projection is None:
        amplitudes=np.zeros(len(e));amplitudes[0]=1.
    else:
        if total_projection<=0:raise ValueError('Projection length must be positive.')
        amplitudes=c*np.exp(-.5*total_projection*(e-e[0]));amplitudes/=np.linalg.norm(amplitudes)
    wave=b@(v@amplitudes);probability=wave*wave
    _,_,px,_,_,diagonals,_,_=basis_data(lx,ly)
    values={key:float(data@probability) for key,data in diagonals.items() if key in KEYS}
    values['ct0']=float((px*px).mean(1)@probability)
    values['mx']=float(np.sign(hx)*(wave@(-transverse@wave)/n))
    for power,key in [(1,'E'),(2,'E2'),(4,'E4')]:values[key]=float(np.sum(amplitudes**2*(e/n)**power))
    values['E_count']=values['E']
    values.update({key:values[original] for key,original in RAW_ALIASES.items()})
    fidelity=float(np.clip(amplitudes[0]**2,0,1))
    return dict(values=values,ground_energy_per_spin=float(e[0]/n),ground_fidelity=fidelity,ground_infidelity=float(np.sum(amplitudes[1:]**2)),
                numerical_residual=residual,dimension=1<<n,block_dimension=len(e),
                projection_total=total_projection,trial='uniform positive product in the |h| gauge',
                reference='unique Perron-Frobenius ground state' if total_projection is None else 'exact normalized two-sided projected pure state',
                sampling_error=0.,sector_gap=float(e[1]-e[0]) if len(e)>1 else None,
                energy_variance=float(np.sum(amplitudes**2*e**2)-np.sum(amplitudes**2*e)**2))


def orbit_probabilities(params,total_projection=None):
    """Probabilities on complete symmetry orbits; no inferred phase/sign data."""
    lx,ly=int(params['Lx']),int(params['Ly'])
    b,_,_,_,e,v,c,_=prepared(lx,ly,float(params['J0']),float(params['J1']),float(params['J2']),abs(float(params['hx'])))
    if total_projection is None:wave=v[:,0]
    else:
        amplitudes=c*np.exp(-.5*total_projection*(e-e[0]));amplitudes/=np.linalg.norm(amplitudes);wave=v@amplitudes
    probability=wave*wave;probability/=probability.sum()
    mapping=b.tocsr().indices.astype('<i4')
    return mapping,probability


def reconstructed_diagnostics(params,probability):
    """A positive, symmetry-invariant state reconstructed from sampled orbits."""
    lx,ly=int(params['Lx']),int(params['Ly']);n=lx*ly
    b,h,transverse,_,e,v,_,residual=prepared(lx,ly,float(params['J0']),float(params['J1']),float(params['J2']),abs(float(params['hx'])))
    psi=np.sqrt(probability/probability.sum());wave=b@psi;energy=float(wave@(h@wave));gap=float(e[1]-e[0])
    return dict(variational_energy_per_spin=energy/n,energy_excess_per_spin=(energy-e[0])/n,
        reconstructed_mx=float(np.sign(params['hx'])*(wave@(-transverse@wave)/n)),
        direct_ground_overlap=float((v[:,0]@psi)**2),
        gap_fidelity_bound_estimate=max(0.,1-max(0.,energy-e[0]+100*float(np.sqrt(len(e)))*residual)/gap),
        bound_note='Variational bound within the invariant sector, using numerical ED gap/residual; not a finite-sample confidence statement.')
