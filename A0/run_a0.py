"""Phase A0 driver: grid gamma x phi x P x T_u of free flames, one JSON per case.

  python run_a0.py run      --mech ../mechanisms/mei2021/mei2021.yaml --tag mei2021 [--no-soret] [-j 56]
  python run_a0.py validate --mech ../mechanisms/mei2021/mei2021.yaml [-j 40]
  python run_a0.py collect  --tag mei2021
  python run_a0.py retry    --mech ... --tag mei2021 [--no-soret]   (failed cases, continuation in phi)

Cases already done (JSON present) are skipped, so a job can simply be resubmitted.
"""
import argparse
import itertools
import json
import os
import signal
import sys
import traceback
from multiprocessing import Pool

import numpy as np

import a0lib

GAMMA = [0.0, 0.1, 0.2, 0.31, 0.43, 0.6, 0.75, 1.0]
PHI = [0.4, 0.5, 0.6, 0.7, 0.8, 1.0]
P_ATM = [1.0, 5.0, 10.0]
T_U = [300.0, 500.0, 750.0]
CASE_TIMEOUT = 3600  # s


def case_name(g, phi, P, T):
    return f"g{g:.2f}_phi{phi:.2f}_P{P:04.1f}_T{T:03.0f}"


class Timeout(Exception):
    pass


def _alarm(signum, frame):
    raise Timeout()


def worker(args):
    mech, outdir, g, phi, P, T, soret = args[:7]
    phi_start = args[7] if len(args) > 7 else None
    name = case_name(g, phi, P, T)
    out = os.path.join(outdir, "cases", name + ".json")
    if os.path.exists(out) and phi_start is None:
        return name, "skip"
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(CASE_TIMEOUT)
    try:
        res = a0lib.run_case(mech, g, phi, P, T, soret=soret, phi_start=phi_start,
                             profile_path=os.path.join(outdir, "profiles", name + ".csv"))
    except Timeout:
        res = {"gamma": g, "phi": phi, "P_atm": P, "T_u": T, "soret": soret,
               "status": "timeout"}
    except Exception as e:  # noqa: BLE001
        res = {"gamma": g, "phi": phi, "P_atm": P, "T_u": T, "soret": soret,
               "status": "failed", "error": str(e)[-500:],
               "trace": traceback.format_exc()[-1500:]}
    finally:
        signal.alarm(0)
    with open(out, "w") as fh:
        fh.write(a0lib.to_json(res))
    return name, res["status"]


def cmd_run(a):
    outdir = os.path.join("results", a.tag + ("" if a.soret else "_nosoret"))
    for d in ("cases", "profiles"):
        os.makedirs(os.path.join(outdir, d), exist_ok=True)
    mech = os.path.abspath(a.mech)
    # hardest (lean, NH3-rich, high P) first for load balance
    grid = sorted(itertools.product(GAMMA, PHI, P_ATM, T_U),
                  key=lambda c: (c[0], c[1], -c[2], c[3]))
    tasks = [(mech, outdir, g, phi, P, T, a.soret) for g, phi, P, T in grid]
    print(f"{len(tasks)} cases -> {outdir}", flush=True)
    with Pool(a.jobs, maxtasksperchild=1) as pool:
        for i, (name, st) in enumerate(pool.imap_unordered(worker, tasks), 1):
            print(f"[{i}/{len(tasks)}] {name} {st}", flush=True)
    cmd_collect(argparse.Namespace(tag=os.path.basename(outdir)))


def cmd_retry(a):
    """Failed or missing cases: continuation in phi from the nearest converged case."""
    global CASE_TIMEOUT
    CASE_TIMEOUT = a.timeout
    import glob
    outdir = os.path.join("results", a.tag + ("" if a.soret else "_nosoret"))
    mech = os.path.abspath(a.mech)
    cases = {c: "missing" for c in itertools.product(GAMMA, PHI, P_ATM, T_U)}
    for f in glob.glob(os.path.join(outdir, "cases", "*.json")):
        d = json.load(open(f))
        cases[(d["gamma"], d["phi"], d["P_atm"], d["T_u"])] = d["status"]
    tasks = []
    for (g, phi, P, T), st in sorted(cases.items()):
        if st == "ok":
            continue
        # nearest converged phi, richer first (continuation works both ways)
        conv = [p for p in PHI if p != phi and cases.get((g, p, P, T)) == "ok"]
        if conv:
            start = min(conv, key=lambda p: (abs(p - phi), p < phi))
            tasks.append((mech, outdir, g, phi, P, T, a.soret, start))
        else:
            print(f"no converged neighbour for {case_name(g, phi, P, T)}")
    print(f"{len(tasks)} cases to retry", flush=True)
    with Pool(max(1, min(a.jobs, len(tasks))), maxtasksperchild=1) as pool:
        for name, st in pool.imap_unordered(worker, tasks):
            print(f"retry {name} {st}", flush=True)
    cmd_collect(argparse.Namespace(tag=os.path.basename(outdir)))


def flatten(d, prefix=""):
    out = {}
    for k, v in d.items():
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            out.update(flatten(v, key + "."))
        else:
            out[key] = v
    return out


def cmd_collect(a):
    import glob
    import pandas as pd
    outdir = os.path.join("results", a.tag)
    rows = [flatten(json.load(open(f))) for f in sorted(glob.glob(os.path.join(outdir, "cases", "*.json")))]
    df = pd.DataFrame(rows).drop(columns=[c for c in ("trace", "mech") if c in rows[0]], errors="ignore")
    df = df.drop(columns=[c for c in df.columns if c.startswith("X_u.")])
    df.to_csv(os.path.join(outdir, "a0_summary.csv"), index=False)
    print(df["status"].value_counts().to_string())
    print("written", os.path.join(outdir, "a0_summary.csv"))


# ---------------------------------------------------------------------------
# validation against Mei et al. 2021 (CNF 231, 111472), supplementary SMM#2
EXP = {  # (gamma, phi, P_atm): (S_L cm/s, err)   T_u = 298 K
    **{(0.4, p, 1.0): v for p, v in zip([0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4],
       [(20.16, .97), (26.26, 1.04), (32.38, 1.13), (35.66, 1.13), (38.11, 1.19), (35.81, 1.13), (30.65, 1.02), (24.82, .93)])},
    **{(0.4, p, 2.0): v for p, v in zip([0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4],
       [(18.54, .93), (22.13, .83), (27.41, .97), (30.27, 1.08), (26.73, .90), (23.48, .85), (19.21, .78)])},
    **{(0.4, p, 5.0): v for p, v in zip([1.1, 1.2, 1.3, 1.4],
       [(19.48, .76), (17.02, .71), (14.78, .68), (12.38, .66)])},
    **{(0.4, p, 10.0): v for p, v in zip([1.1, 1.2, 1.3, 1.4],
       [(13.66, 1.21), (10.98, .62), (9.70, .61), (7.96, .62)])},
    **{(g, 1.0, 1.0): v for g, v in zip([0.1, 0.2, 0.3, 0.5, 0.6, 0.7, 0.8],
       [(9.43, .83), (14.71, .81), (24.20, .90), (56.60, 2.67), (68.27, 2.06), (90.06, 2.48), (112.66, 3.06)])},
    **{(g, 1.0, 2.0): v for g, v in zip([0.1, 0.2, 0.3, 0.5, 0.6],
       [(7.09, .73), (11.36, .72), (17.62, .82), (41.37, 1.50), (57.90, 1.84)])},
    **{(g, 1.0, 5.0): v for g, v in zip([0.1, 0.2, 0.3],
       [(5.70, .65), (7.54, .65), (12.03, .68)])},
}


def _validate_one(args):
    mech, (g, phi, P), soret = args
    import cantera as ct
    gas = ct.Solution(mech)
    a0lib.set_mixture(gas, g, phi, 298.0, P)
    f = a0lib.solve_free_flame(gas, a0lib.initial_width(P, 298.0), soret=soret)
    return g, phi, P, soret, f.velocity[0] * 100


def cmd_validate(a):
    import pandas as pd
    mech = os.path.abspath(a.mech)
    tasks = [(mech, k, s) for k in EXP for s in (False, True)]
    with Pool(a.jobs) as pool:
        out = pool.map(_validate_one, tasks)
    df = pd.DataFrame(out, columns=["gamma", "phi", "P_atm", "soret", "sL_sim"])
    df["sL_exp"] = [EXP[(g, p, P)][0] for g, p, P in zip(df.gamma, df.phi, df.P_atm)]
    df["err_exp"] = [EXP[(g, p, P)][1] for g, p, P in zip(df.gamma, df.phi, df.P_atm)]
    df["rel_dev"] = df.sL_sim / df.sL_exp - 1
    os.makedirs("results", exist_ok=True)
    tag = os.path.splitext(os.path.basename(mech))[0]
    df.to_csv(f"results/validation_{tag}.csv", index=False)
    print(df.to_string(index=False, float_format="%.3f"))
    for s, d in df.groupby("soret"):
        print(f"soret={s}: mean |dev| = {np.mean(np.abs(d.rel_dev)):.3f}, max = {np.max(np.abs(d.rel_dev)):.3f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run", "validate", "collect", "retry"])
    ap.add_argument("--mech", default="../mechanisms/mei2021/mei2021.yaml")
    ap.add_argument("--tag", default="mei2021")
    ap.add_argument("--no-soret", dest="soret", action="store_false")
    ap.add_argument("-j", "--jobs", type=int, default=int(os.environ.get("SLURM_CPUS_ON_NODE", 4)))
    ap.add_argument("--timeout", type=int, default=600, help="retry: s per case")
    a = ap.parse_args()
    {"run": cmd_run, "validate": cmd_validate, "collect": cmd_collect, "retry": cmd_retry}[a.cmd](a)
    sys.exit(0)
