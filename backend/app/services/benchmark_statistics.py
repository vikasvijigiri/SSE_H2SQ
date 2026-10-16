# Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
"""Chain-aware FFT autocorrelations, blocked jackknife and rank-split R-hat."""
import math
import numpy as np
from scipy.stats import rankdata, norm

PRIMARY = ["E", "E2", "E4", "mz2", "mz4", "mx", "plaquette", "ice", "rho", "stripe2", "polar2", "ct0", "ct1", "ct2", "ct3", "ct4"]


def autocorrelation(samples):
    x = np.asarray(samples, dtype=float)
    x = x-x.mean()
    if np.max(np.abs(x)) < 1e-14:
        return {"tau": None, "tau_error": None, "ess": None, "acf": None, "window": 0, "constant": True}
    n = len(x)
    size = 1 << (2*n-1).bit_length()
    spectrum = np.fft.rfft(x, size)
    covariance = np.fft.irfft(spectrum*spectrum.conj(), size)[:n] / np.arange(n, 0, -1)
    acf = covariance/covariance[0]
    # Initial positive, monotone paired sequence. No chain boundaries are joined.
    total, previous, window = 0., float("inf"), 0
    for lag in range(1, min(n-1, n//5, 10000), 2):
        pair = min(previous, float(acf[lag]+acf[lag+1]))
        if pair <= 0:
            break
        total += pair
        previous, window = pair, lag+1
    tau = max(.5, .5+total)
    return {"tau": tau, "tau_error": tau*math.sqrt((4*window+2)/n), "ess": n/(2*tau), "acf": acf, "window": window, "constant": False}


def split_rhat(chains):
    half = min(map(len, chains))//2
    split = np.asarray([part for chain in chains for part in (chain[:half], chain[-half:])])
    def rhat(values):
        within = np.var(values, axis=1, ddof=1).mean()
        between = half*np.var(values.mean(axis=1), ddof=1)
        if within < 1e-28:
            return None if between < 1e-28 else float("inf")
        return float(np.sqrt(((half-1)*within/half+between/half)/within))
    if np.ptp(split) < 1e-14:
        return None
    def normalize(values):
        ranks = rankdata(values.ravel(), method="average")
        return norm.ppf((ranks-.375)/(len(ranks)+.25)).reshape(values.shape)
    values = [rhat(normalize(split)), rhat(normalize(np.abs(split-np.median(split))))]
    return max(value for value in values if value is not None)


def derived(means, beta, n, molecules, hx):
    values = {key: float(means[key]) for key in PRIMARY}
    values["cv"] = beta**2*n*(means["E2"]-means["E"]**2)
    values["chi_stag"] = beta*n*(means["bar_m2"]-means["bar_m"]**2)
    values["chi_polar"] = beta*molecules*(means["bar_p2"]-means["bar_px"]**2-means["bar_py"]**2)
    values["chi_x"] = (means["nxf2"]-means["nx"]**2)/(beta*hx**2*n) if hx else None
    values["binder"] = .5*(1-means["stripe4"]/(3*means["stripe2"]**2)) if means["stripe2"] > 1e-16 else None
    return values


def influence(chains, means, key, beta, n, molecules, hx):
    if key in PRIMARY:
        return [chain[key] for chain in chains]
    gradients = {
        "cv": {"E2": beta**2*n, "E": -2*beta**2*n*means["E"]},
        "chi_stag": {"bar_m2": beta*n, "bar_m": -2*beta*n*means["bar_m"]},
        "chi_polar": {"bar_p2": beta*molecules, "bar_px": -2*beta*molecules*means["bar_px"], "bar_py": -2*beta*molecules*means["bar_py"]},
    }
    if hx:
        gradients["chi_x"] = {"nxf2": 1/(beta*hx**2*n), "nx": -2*means["nx"]/(beta*hx**2*n)}
    if means["stripe2"] > 1e-16:
        gradients["binder"] = {"stripe4": -1/(6*means["stripe2"]**2), "stripe2": means["stripe4"]/(3*means["stripe2"]**3)}
    return [sum(coefficient*chain[column] for column, coefficient in gradients[key].items()) for chain in chains]


def summarize(chains, params):
    beta, n = params["beta"], params["Lx"]*params["Ly"]
    molecules, hx = n//2, params["hx"]
    columns = chains[0].dtype.names
    sizes=np.asarray([len(chain) for chain in chains],dtype=float)
    sample_count=int(sizes.sum())
    chain_weights=sizes/sample_count
    chain_means=[{key:float(chain[key].mean()) for key in columns} for chain in chains]
    full_means={key:float(sum(weight*means[key] for weight,means in zip(chain_weights,chain_means))) for key in columns}
    chain_estimates=[derived(means,beta,n,molecules,hx) for means in chain_means]
    values = derived(full_means, beta, n, molecules, hx)
    diagnostics, projected = {}, {}
    largest_tau = .5
    for key, value in values.items():
        if value is None:
            continue
        streams = influence(chains, full_means, key, beta, n, molecules, hx)
        projected[key] = streams
        stats = [autocorrelation(stream) for stream in streams]
        taus = [stat["tau"] for stat in stats if stat["tau"] is not None]
        if taus:
            largest_tau = max(largest_tau, max(taus))
        diagnostics[key] = stats
    length = min(map(len, chains))
    requested_block = max(32, math.ceil(10*largest_tau))
    block_size = min(requested_block, max(1, length//8))
    block_counts=[len(chain)//block_size for chain in chains]
    blocks = np.concatenate([np.column_stack([chain[key][:count*block_size].reshape(count,block_size).mean(axis=1) for key in columns]) for chain,count in zip(chains,block_counts)])
    # All rows enter the point estimate, including the final partial block.
    means=full_means
    values = derived(means, beta, n, molecules, hx)
    leave_out = []
    sample_sum=np.asarray([means[key]*sample_count for key in columns])
    for block in blocks:
        leave_out.append(derived(dict(zip(columns,(sample_sum-block_size*block)/(sample_count-block_size))),beta,n,molecules,hx))
    result = {}
    for key, value in values.items():
        if value is None:
            result[key] = {"mean": None, "error": None, "status": "estimator unavailable"}
            continue
        stats = diagnostics[key]
        streams = projected[key]
        rhat = split_rhat(streams) if len(chains) > 1 else None
        estimates = [row[key] for row in leave_out]
        valid = all(estimate is not None and math.isfinite(estimate) for estimate in estimates)
        # Equal-size delete-block jackknife, normalized to the actual row
        # count when final partial blocks are held in the point estimate.
        jackknife_error = math.sqrt((sample_count-block_size)**2/(sample_count*block_size)*float(np.var(estimates,ddof=1))) if valid else None
        chain_values = [estimate[key] for estimate in chain_estimates]
        center=float(np.dot(chain_weights,chain_values))
        between_error=math.sqrt(float(np.sum(chain_weights**2*(np.asarray(chain_values)-center)**2/(1-chain_weights)))) if len(chains)>1 else 0.
        gamma_variance = sum(len(stream)*np.var(stream,ddof=1)*2*stat["tau"] for stream,stat in zip(streams,stats) if stat["tau"] is not None)/sample_count**2
        error = max(jackknife_error or 0., between_error, math.sqrt(gamma_variance))
        constant = all(stat["constant"] for stat in stats)
        analytic_zero = key == "mx" and hx == 0
        flags = []
        if rhat is not None and (not math.isfinite(rhat) or rhat > 1.01):
            flags.append("between-chain disagreement")
        if constant and not analytic_zero:
            flags.append("constant series: mixing unresolved")
            error = max(error, between_error) if between_error > 0 else None
        ess = sum(stat["ess"] or 0 for stat in stats)
        if not constant and ess < 400:
            flags.append("low effective sample size")
        if requested_block > length//8:
            flags.append("chain too short for correlation scale")
        if len(chains) < 4:
            flags.append("fewer than four independent chains")
        lag_indices = np.unique(np.concatenate(([0], np.geomspace(1, min(length//5, 2048), 64).astype(int))))
        acfs = [stat["acf"] for stat in stats if stat["acf"] is not None]
        curve = [{"lag": int(lag), "value": float(np.mean([acf[lag] for acf in acfs]))} for lag in lag_indices] if acfs else []
        taus = [stat["tau"] for stat in stats if stat["tau"] is not None]
        result[key] = {"mean": float(value), "error": float(error) if error is not None else None, "tau_int": float(np.mean(taus)) if taus else None, "tau_error": float(np.mean([stat["tau_error"] for stat in stats if stat["tau_error"] is not None])) if taus else None, "ess": float(ess) if taus else None, "rhat": rhat if rhat is not None and math.isfinite(rhat) else None, "acf": curve, "flags": flags, "status": "analytic zero" if analytic_zero else "review" if flags else "checks passed"}
    comparisons={}
    for key in ("mz2","stripe2","polar2","ct0","ct2","ct4"):
        raw="boundary_"+key
        if raw not in columns:continue
        streams=[chain[raw] for chain in chains]
        raw_stats=[autocorrelation(stream) for stream in streams]
        mean=full_means[raw]
        raw_blocks=np.concatenate([stream[:count*block_size].reshape(count,block_size).mean(axis=1) for stream,count in zip(streams,block_counts)])
        block_error=float(np.std(raw_blocks,ddof=1)*np.sqrt(block_size/sample_count))
        acf_error=math.sqrt(sum(len(stream)*np.var(stream,ddof=1)*2*stat["tau"] for stream,stat in zip(streams,raw_stats) if stat["tau"] is not None)/sample_count**2)
        raw_chain_means=np.asarray([values[raw] for values in chain_means])
        chain_error=math.sqrt(float(np.sum(chain_weights**2*(raw_chain_means-mean)**2/(1-chain_weights)))) if len(chains)>1 else 0.
        raw_error=max(block_error,acf_error,chain_error)
        if all(stat["constant"] for stat in raw_stats) and raw_error==0:raw_error=None
        improved=result[key]
        comparisons[key]={"boundary_mean":mean,"boundary_SE":raw_error,"improved_mean":improved["mean"],"improved_SE":improved["error"],"SE_gain":raw_error/improved["error"] if raw_error and improved["error"] else None}
    return {"observables": result,"estimator_comparison":comparisons,"block_size":block_size,"blocks":len(blocks),"retained_samples":sample_count,"samples_in_full_blocks":sum(block_counts)*block_size,"partial_block_rows_per_chain":[len(chain)-count*block_size for chain,count in zip(chains,block_counts)],"discarded_tail_per_chain":0,"chains":len(chains),"error_method":"sample-count weighted means; maximum of delete-block jackknife, chain-wise ACF variance and HC2 independent-chain SE; every row enters the point estimate","tau_convention":"tau_int=1/2+sum rho(t), in complete sweeps; derived observables use linearized influence series"}
