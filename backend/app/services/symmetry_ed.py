# Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
"""Complete 4x4 eigenspectrum in commuting real symmetry sectors.

Half translations, reflections through a plaquette centre and spin inversion
generate an abelian group of involutions. No magnetization-sector restriction
is made: the transverse field mixes magnetization. All characters and all
states in every allowed character sector contribute to the thermal trace.
"""
from functools import lru_cache
from pathlib import Path
import hashlib
import inspect
import json
import os
import time

import numpy as np
from scipy import sparse
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / ".cache" / "symmetry-ed"
VERSION = "complete-real-reflection-sectors-v1"


def _permutation(indices, mapping):
    result = np.zeros(len(indices), dtype=np.int32)
    for source, target in enumerate(mapping):
        result |= ((indices >> source) & 1).astype(np.int32) << target
    return result


@lru_cache(maxsize=2)
def sectors(lx=4, ly=4):
    if (lx, ly) not in {(2, 2), (4, 4)}:
        raise ValueError("Real reflection sectors support 2x2 tests and 4x4 production.")
    n, dimension = lx*ly, 1 << (lx*ly)
    indices = np.arange(dimension, dtype=np.int32)
    maps = []
    for transform in (lambda x,y: ((x+lx//2) % lx,y),
                      lambda x,y: (x,(y+ly//2) % ly),
                      lambda x,y: ((1-x) % lx,y),
                      lambda x,y: (x,(1-y) % ly)):
        maps.append(_permutation(indices, [tx+lx*ty for x,y in
                    ((site % lx, site//lx) for site in range(n))
                    for tx,ty in [transform(x,y)]]))
    maps.append(indices ^ (dimension-1))
    for a in maps:
        if not np.array_equal(a[a], indices):
            raise RuntimeError("Symmetry generator is not an involution.")
        for b in maps:
            if not np.array_equal(a[b], b[a]):
                raise RuntimeError("Symmetry generators do not commute.")
    group = np.empty((32, dimension), dtype=np.int32)
    group[0] = indices
    for element in range(1,32):
        bit = (element & -element).bit_length()-1
        group[element] = maps[bit][group[element ^ (1 << bit)]]
    representatives = group.min(axis=0)
    reps = np.unique(representatives)
    characters = np.array([[(-1)**((character & element).bit_count() % 2)
                            for element in range(32)] for character in range(32)])
    result = []
    for character in range(32):
        rows, cols, data, orbit_ids = [], [], [], []
        for orbit_id, representative in enumerate(reps):
            images = group[:, representative]
            # Stabilizer characters must be trivial for a nonzero projection.
            if np.any(characters[character, images == representative] != 1):
                continue
            states, first = np.unique(images, return_index=True)
            column = len(orbit_ids)
            rows.extend(states); cols.extend([column]*len(states))
            data.extend(characters[character, first]/np.sqrt(len(states)))
            orbit_ids.append(orbit_id)
        basis = sparse.csc_matrix((data,(rows,cols)),shape=(dimension,len(orbit_ids)))
        if basis.shape[1]:
            result.append({"character": character,"basis":basis,
                           "orbits":np.asarray(orbit_ids,dtype=np.int32)})
    if sum(item["basis"].shape[1] for item in result) != dimension:
        raise RuntimeError("Symmetry sectors do not span the full Hilbert space.")
    return result, reps, group


def hamiltonian(diagonal, hx, n):
    dimension = len(diagonal)
    indices = np.arange(dimension)
    rows = np.repeat(indices,n)
    cols = (indices[:,None] ^ (1 << np.arange(n))).ravel()
    transverse = sparse.csr_matrix((np.full(len(rows),-1.),(rows,cols)),shape=(dimension,dimension))
    return sparse.diags(diagonal,format="csr") + hx*transverse, transverse


def build_spectrum(lx,ly,j0,j1,j2,hx,diagonal,observables):
    """Persist exact sector eigenpairs so beta sweeps reuse one diagonalization."""
    fingerprint = hashlib.sha256((VERSION+str((lx,ly,j0,j1,j2,hx))+
                 inspect.getsource(sectors.__wrapped__)+inspect.getsource(_permutation)+
                 inspect.getsource(hamiltonian)).encode()+diagonal.tobytes()).hexdigest()[:24]
    folder = CACHE / fingerprint
    folder.mkdir(parents=True,exist_ok=True)
    blocks,reps,group = sectors(lx,ly)
    h, transverse = hamiltonian(diagonal,hx,lx*ly)
    # Test the supplied Hamiltonian, not just a presumed lattice symmetry.
    for permutation in group[[1,2,4,8,16]]:
        if not np.allclose(diagonal, diagonal[permutation],atol=1e-12,rtol=0):
            raise RuntimeError("Hamiltonian does not preserve the selected spatial symmetries.")
    residual = 0.
    output = []
    for sector in blocks:
        character,basis = sector["character"],sector["basis"]
        path=folder/f"sector-{character}.npz"
        if path.exists():
            with np.load(path,allow_pickle=False) as saved:
                energies=saved["energies"];vectors=saved["vectors"]
                local_residual=float(saved["residual"])
        else:
            started=time.monotonic()
            block=(basis.T @ h @ basis).toarray()
            if np.max(np.abs(block-block.T)) > 1e-11:
                raise RuntimeError("Projected Hamiltonian is not Hermitian.")
            energies,vectors=eigh(block,driver="evd",check_finite=False,overwrite_a=False)
            local_residual=float(np.max(np.linalg.norm(block @ vectors-vectors*energies,axis=0)))
            temporary=folder/f"sector-{character}.tmp.npz"
            np.savez(temporary,energies=energies,vectors=vectors,residual=local_residual)
            temporary.replace(path)
            print(f"ED sector {character}: {len(energies)} states, {time.monotonic()-started:.1f}s",flush=True)
        residual=max(residual,local_residual)
        # Observable vectors are diagonal in the original spin basis.
        diagonal_means=np.asarray(basis.power(2).T @ np.column_stack(list(observables.values())))
        expectations=diagonal_means.T @ (vectors*vectors)
        projected_x=basis.T @ (-transverse/(lx*ly)) @ basis
        mx=np.sum(vectors*(projected_x @ vectors),axis=0)
        output.append({**sector,"energies":energies,"path":path,
                       "expectations":dict(zip(observables,expectations)),"mx":mx})
        del vectors
    if sum(len(item["energies"]) for item in output) != 1 << (lx*ly):
        raise RuntimeError("Incomplete eigenspectrum; refusing a thermal reference.")
    return output,residual,folder


def _components(operator, group):
    """Character-project a diagonal operator, retaining nonzero components."""
    characters=np.array([[(-1)**((a&b).bit_count() % 2) for b in range(32)] for a in range(32)])
    projected=characters @ operator[group]/32
    return {q: values for q,values in enumerate(projected) if np.max(np.abs(values)) > 1e-12}


def _vectors(sector):
    with np.load(sector["path"],allow_pickle=False) as saved:
        return saved["vectors"]


def response_reference(blocks,beta,ground,partition,diagonal_ops,transverse,folder,
                       means,n,molecules,tail_tolerance=1e-13):
    """Lehmann responses with all intermediate states and bounded thermal tail.

    The complete spectrum is used for Z and every equal-time expectation.
    Response matrix elements are omitted only if BOTH endpoint states belong
    to a thermal tail whose exact total probability is <= tail_tolerance.
    Convexity bounds the missing positive c(tau) by p_tail ||O||^2, and a Kubo
    integral by beta times that bound. High-energy intermediate states linked
    to a retained endpoint are never dropped.
    """
    analysis_hash=hashlib.sha256((inspect.getsource(response_reference)+inspect.getsource(_components)).encode()).hexdigest()[:16]
    path=folder/f"response-{beta:.17g}-{analysis_hash}.json"
    if path.exists():
        return json.loads(path.read_text())
    energies=np.sort(np.concatenate([block["energies"] for block in blocks]))
    probabilities=np.exp(-beta*(energies-ground))/partition
    reverse=np.cumsum(probabilities[::-1])[::-1]
    eligible=np.flatnonzero(reverse <= tail_tolerance)
    cutoff=energies[eligible[0]-1] if len(eligible) and eligible[0] else energies[-1]
    # Include an entire degenerate eigenspace at the boundary.
    active={b["character"]:b["energies"] <= cutoff+1e-10 for b in blocks}
    tail=sum(float(np.exp(-beta*(b["energies"][~active[b['character']]]-ground)).sum()/partition) for b in blocks)
    _,_,group=sectors(int(np.sqrt(n)),int(np.sqrt(n)))
    channels={}
    norms={}
    for name,operator,normalization,kind in diagonal_ops:
        norms[name]=norms.get(name,0.)+normalization*float(np.max(np.abs(operator)))**2
        for q,component in _components(operator,group).items():
            channels.setdefault(q,[]).append((name,component,normalization,kind))
    channels.setdefault(0,[]).append(("chi_x",None,n,"kubo"))
    norms["chi_x"]=n
    output={"chi_stag":0.,"chi_polar":0.,"chi_x":0.,**{f"ct{i}":0. for i in range(1,5)}}
    lookup={b["character"]:b for b in blocks}
    fractions=[.125,.25,.375,.5]
    for left in blocks:
        a=left["character"]
        va=_vectors(left)
        ia=np.flatnonzero(active[a]);inactive_a=np.flatnonzero(~active[a])
        for q,operators in channels.items():
            b=a^q
            if b<a or b not in lookup:
                continue
            right=lookup[b];vb=va if a==b else _vectors(right)
            ib=np.flatnonzero(active[b]);inactive_b=np.flatnonzero(~active[b])
            if not len(ia) and not len(ib):
                continue
            for name,operator,normalization,kind in operators:
                basis_a,basis_b=left["basis"],right["basis"]
                matrix=(basis_a.T @ (transverse @ basis_b) if operator is None else
                        basis_a.T @ basis_b.multiply(operator[:,None]))
                # Each unordered cross-sector pair is covered exactly once.
                rectangles=[]
                if len(ia):
                    elements=va[:,ia].T @ (matrix @ vb)
                    rectangles.append((ia,np.arange(len(right["energies"])),elements,a==b))
                if a!=b and len(inactive_a) and len(ib):
                    elements=(vb[:,ib].T @ (matrix.T @ va[:,inactive_a])).T
                    rectangles.append((inactive_a,ib,elements,False))
                for rows,cols,elements,same in rectangles:
                    ea=left["energies"][rows]-ground;eb=right["energies"][cols]-ground
                    amplitude=normalization*(elements*elements)
                    extra=(~active[b][cols]).astype(float)[None,:] if same else 1.
                    if kind=="kubo":
                        pa=np.exp(-beta*ea)[:,None]/partition;pb=np.exp(-beta*eb)[None,:]/partition
                        delta=eb[None,:]-ea[:,None]
                        distance=np.abs(delta)
                        maximum=np.maximum(pa,pb)
                        kernel=np.divide(maximum*(-np.expm1(-beta*distance)),distance,
                                         out=beta*maximum.copy(),where=distance>1e-13)
                        output[name]+=float(np.sum(amplitude*kernel*(1+extra)))
                    else:
                        for i,fraction in enumerate(fractions,1):
                            forward=np.exp(-beta*((1-fraction)*ea[:,None]+fraction*eb[None,:]))/partition
                            reverse_weight=np.exp(-beta*(fraction*ea[:,None]+(1-fraction)*eb[None,:]))/partition
                            output[f"ct{i}"]+=float(np.sum(amplitude*(forward+extra*reverse_weight)))
            if a!=b:
                del vb
        del va
    output["chi_x"]-=beta*n*means["mx"]**2
    bounds={name:beta*tail*norms[name] for name in ("chi_stag","chi_polar","chi_x")}
    bounds.update({f"ct{i}":tail*norms["ct"] for i in range(1,5)})
    result={"values":output,"absolute_error_bounds":bounds,"omitted_thermal_probability":tail,
            "retained_thermal_states":sum(int(v.sum()) for v in active.values()),
            "intermediate_states":"all 65536 eigenstates; only pairs with both endpoints in the bounded thermal tail are omitted"}
    temporary=path.with_suffix(".tmp")
    temporary.write_text(json.dumps(result,allow_nan=False));temporary.replace(path)
    return result


@lru_cache(maxsize=2)
def _prepared(lx,ly,j0,j1,j2,hx):
    from .physics_model import basis_data
    _,local,px,py,pairs,observables,polx,poly=basis_data(lx,ly)
    diagonal=-j0*np.prod(local,axis=2).sum(axis=1)+j1*(local[:,:,0]*local[:,:,2]+local[:,:,1]*local[:,:,3]).sum(axis=1)
    for a,b in pairs:
        diagonal=diagonal+j2*(px[:,a]*px[:,b]+py[:,a]*py[:,b])
    all_observables={**observables,"ct0":(px*px).mean(axis=1)}
    blocks,residual,folder=build_spectrum(lx,ly,j0,j1,j2,hx,diagonal,all_observables)
    return blocks,residual,folder,diagonal,px,observables,polx,poly


def exact_4x4_reference(params):
    if float(params["hx"])<0:
        # The product of all sigma_z operators changes h -> -h and leaves
        # every diagonal interaction invariant. Reuse the identical spectrum.
        result=json.loads(json.dumps(exact_4x4_reference({**params,"hx":-float(params["hx"])})))
        result["values"]["mx"]*=-1
        return result
    lx,ly=int(params["Lx"]),int(params["Ly"])
    n=lx*ly;beta=float(params["beta"])
    blocks,residual,folder,diagonal,px,observables,polx,poly=_prepared(lx,ly,float(params["J0"]),float(params["J1"]),float(params["J2"]),float(params["hx"]))
    ground=min(float(b["energies"].min()) for b in blocks)
    partition=sum(float(np.exp(-beta*(b["energies"]-ground)).sum()) for b in blocks)
    values={key:0. for key in blocks[0]["expectations"]}
    values.update(E=0.,E2=0.,E4=0.,mx=0.)
    for block in blocks:
        weight=np.exp(-beta*(block["energies"]-ground))/partition
        for name,expectation in block["expectations"].items():
            values[name]+=float(expectation @ weight)
        for key,power in (("E",1),("E2",2),("E4",4)):
            values[key]+=float((block["energies"]/n)**power @ weight)
        values["mx"]+=float(block["mx"] @ weight)
    values["cv"]=max(0.,beta*beta*n*(values["E2"]-values["E"]**2))
    values["binder"]=.5*(1-values["stripe4"]/(3*values["stripe2"]**2)) if values["stripe2"]>0 else None
    _,transverse=hamiltonian(diagonal,0,n)
    # The two plaquette sublattices are symmetry-related by a (1,1) translation.
    # Keep both representatives so this identity is not assumed in the response.
    other=2 if n==16 else 1
    operators=[("chi_stag",observables["mz"],n,"kubo"),
               ("chi_polar",polx,n//2,"kubo"),("chi_polar",poly,n//2,"kubo"),
               ("ct",px[:,0],.5,"correlation"),("ct",px[:,other],.5,"correlation")]
    response=response_reference(blocks,beta,ground,partition,operators,-transverse/n,folder,values,n,n//2)
    values.update(response["values"])
    return {"values":values,"statistical_errors":{key:0. for key in values},
            "numerical_residual":residual,"absolute_error_bounds":response["absolute_error_bounds"],
            "reference":"complete symmetry-resolved eigenspectrum; bounded Lehmann thermal tail",
            "dimension":sum(len(b["energies"]) for b in blocks),"sector_dimensions":[len(b["energies"]) for b in blocks],
            "ground_energy_per_spin":ground/n,"thermal_response":response,
            "model_version":"cyclic-pair-j1-sublattice-dipoles-v3"}
