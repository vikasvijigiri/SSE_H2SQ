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
    full_means = {key: float(np.mean([np.mean(chain[key]) for chain in chains])) for key in columns}
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
    block_count = length//block_size
    blocks = np.asarray([[np.mean(chain[key][:block_count*block_size].reshape(block_count, block_size), axis=1) for key in columns] for chain in chains])
    blocks = blocks.transpose(0, 2, 1).reshape(-1, len(columns))
    means = dict(zip(columns, blocks.mean(axis=0)))
    values = derived(means, beta, n, molecules, hx)
    leave_out = []
    block_sum = blocks.sum(axis=0)
    for block in blocks:
        leave_out.append(derived(dict(zip(columns, (block_sum-block)/(len(blocks)-1))), beta, n, molecules, hx))
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
        jackknife_error = math.sqrt((len(blocks)-1)*float(np.var(estimates))) if valid else None
        chain_values = [derived({column: float(chain[column].mean()) for column in columns}, beta, n, molecules, hx)[key] for chain in chains]
        between_error = float(np.std(chain_values, ddof=1)/np.sqrt(len(chains))) if len(chains) > 1 else 0.
        gamma_variance = sum(np.var(stream, ddof=1)*2*stat["tau"]/len(stream) for stream, stat in zip(streams, stats) if stat["tau"] is not None)/len(chains)**2
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
    return {"observables": result, "block_size": block_size, "blocks": len(blocks), "retained_samples": block_count*block_size*len(chains), "discarded_tail_per_chain": length-block_count*block_size, "chains": len(chains), "error_method": "maximum of blocked jackknife, per-chain ACF variance estimate and independent-chain SE", "tau_convention": "tau_int=1/2+sum rho(t), in complete sweeps; derived observables use linearized influence series"}
