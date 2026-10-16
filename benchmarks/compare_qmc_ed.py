"""
Run the three validation test cases and print a comparison table.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import requests, json

API = "http://127.0.0.1:8001/api"

TESTS = [
    {"name": "J0+J1 only",   "J0": 1.0, "J1": 0.3, "J2": 0.0,   "J3": 0.0, "hx": 0.0},
    {"name": "J0+J1+hx",     "J0": 1.0, "J1": 0.3, "J2": 0.0,   "J3": 0.0, "hx": 0.3},
    {"name": "J0+J1+J2+hx",  "J0": 1.0, "J1": 0.3, "J2": 0.025, "J3": 0.0, "hx": 0.3},
]

COMMON    = {"Lx": 4, "Ly": 4, "beta": 1.0}
QMC_EXTRA = {"Lz": 1, "n_therm": 200000, "n_measure": 200000}
ED_EXTRA  = {"Nl": 1, "n": 50, "sparse": True}

OBS = [
    ("enrg",        "<E>",   "float"),
    ("enrg2",       "<E2>",  "float"),
    ("enrg4",       "<E4>",  "float"),
    ("SMag_square", "<m2>",  "sci"),
    ("SMag_four",   "<m4>",  "sci"),
    ("Mag_x",       "<Mx>",  "sci"),
]

SEP  = "=" * 88
DASH = "-" * 88
THIN = "-" * 88

def fmt(v, kind):
    if v is None:
        return "        --    "
    if kind == "float":
        return f"{v:+.6f}"
    return f"{v:+.4e}"

def run(endpoint, payload):
    try:
        r = requests.post(f"{API}/{endpoint}", json=payload, timeout=1200)
        r.raise_for_status()
        return r.json()
    except Exception as e:
        return {"error": str(e)}

def extract(res, key):
    if not res or "data" not in res or not res["data"]:
        return None
    d = res["data"]
    if isinstance(d, list) and d:
        item = d[0]
        if isinstance(item, dict):
            return item.get(key)
    return None

print()
print(SEP)
print("  SSE-H2SQ Validation: QMC (loop-cluster update) vs Exact Diagonalization")
print(SEP)

for t in TESTS:
    print()
    print(DASH)
    print(f"  Test: {t['name']}  |  J0={t['J0']} J1={t['J1']} J2={t['J2']} hx={t['hx']}  beta=1.0  4x4")
    print(DASH)

    qmc_p = {**COMMON, **QMC_EXTRA, "J0": t["J0"], "J1": t["J1"],
             "J2": t["J2"], "J3": t["J3"], "hx": t["hx"]}
    ed_p  = {**COMMON, **ED_EXTRA,  "J0": t["J0"], "J1": t["J1"],
             "J2": t["J2"], "h": t["hx"]}

    print("  [QMC] running ... ", end="", flush=True)
    qmc_res = run("run-qmc", qmc_p)
    print("done")

    print("  [ED]  running ... ", end="", flush=True)
    ed_res  = run("run-ed",  ed_p)
    print("done")

    if "error" in qmc_res:
        print(f"  QMC ERROR: {qmc_res['error']}")
    if "error" in ed_res:
        print(f"  ED  ERROR: {ed_res['error']}")

    qmc = extract(qmc_res, "qmc")
    ed  = extract(ed_res,  "ed")

    print()
    print(f"  {'Observable':<14} {'ED':>14} {'QMC':>14} {'|delta|':>13}  {'Agree%':>8}")
    print(f"  {'-'*14} {'-'*14} {'-'*14} {'-'*13}  {'-'*8}")
    for key, label, kind in OBS:
        ev = ed.get(key)  if ed  else None
        qv = qmc.get(key) if qmc else None
        ds = ag = ""
        if ev is not None and qv is not None:
            d  = abs(qv - ev)
            ds = f"{d:.4e}"
            ag = f"{max(0, 1 - d/abs(ev))*100:.3f}%" if abs(ev) > 1e-12 else "  N/A"
        print(f"  {label:<14} {fmt(ev,kind):>14} {fmt(qv,kind):>14} {ds:>13}  {ag:>8}")

print()
print(SEP)
print()
