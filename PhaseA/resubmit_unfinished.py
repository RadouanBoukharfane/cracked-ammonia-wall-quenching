"""Finds Phase A runs whose job has ended before amr.stop_time (crash or wall time), archives their logs in
crashN/ and resubmits run.sbatch (restart from the last chk). Updates cases/<--jobs>.

  python resubmit_unfinished.py [--dry] [--glob 'runs/A_*'] [--nice 1000] [--jobs wave2_jobs.txt]
"""
import argparse
import glob
import json
import os
import re
import shutil
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
CLEAN_ENV = "env -i HOME=$HOME USER=$USER LOGNAME=$USER TERM=dumb PATH=/usr/bin:/bin bash -lc"


def active_jobs():
    out = subprocess.run(["squeue", "-u", os.environ["USER"], "-h", "-o", "%i %Z"], capture_output=True, text=True)
    return {line.split()[1] for line in out.stdout.splitlines() if len(line.split()) == 2}


def last_time(run):
    log = os.path.join(run, "pele.log")
    if not os.path.exists(log):
        return None
    t = re.findall(r"STEP \[\d+\] - Time: ([0-9.eE+-]+),", open(log, errors="replace").read())
    return float(t[-1]) if t else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--glob", default="runs/A_*")
    ap.add_argument("--nice", type=int, default=0, help="sbatch --nice value (wave 2: 1000)")
    ap.add_argument("--jobs", default="wave1_jobs.txt", help="job list in cases/ to update")
    a = ap.parse_args()
    busy = active_jobs()
    jobs_file = os.path.join(HERE, "cases", a.jobs)
    for run in sorted(glob.glob(os.path.join(HERE, a.glob))):
        if os.path.abspath(run) in busy or not os.path.exists(os.path.join(run, "pele.log")):
            continue
        m = json.load(open(os.path.join(run, "case.json")))
        stop = float(re.search(r"amr.stop_time = (\S+)", open(os.path.join(run, "input.2d")).read()).group(1))
        t = last_time(run)
        if t is None or t >= 0.999 * stop:
            continue
        # linear window already closed (max A > 50 A0 before the stop): omega(k) is complete, no restart
        f = os.path.join(run, "front.npz")
        if os.path.exists(f):
            import numpy as np
            fr = np.load(f)
            yf = fr["yf"]
            A = 2 * np.abs(np.fft.rfft(yf - yf.mean(1, keepdims=True), axis=1))[:, 1:m["nmodes"] + 1] / yf.shape[1]
            if (A.max(axis=1) > 50 * m["A0"]).any():
                print(f"{os.path.basename(run)}: stopped at {t / m['tau_F']:.2f} tau_F but linear window closed: skip")
                continue
        chk = sorted(glob.glob(os.path.join(run, "chk*")))
        print(f"{os.path.basename(run)}: stopped at {t / m['tau_F']:.2f} tau_F, restart from "
              f"{os.path.basename(chk[-1]) if chk else 'scratch'}")
        if a.dry:
            continue
        n = 1
        while os.path.exists(os.path.join(run, f"crash{n}")):
            n += 1
        arch = os.path.join(run, f"crash{n}")
        os.makedirs(arch)
        for f in glob.glob(os.path.join(run, "Backtrace.*")) + glob.glob(os.path.join(run, "slurm-*.out")) + \
                [os.path.join(run, "pele.log")]:
            shutil.move(f, arch)
        j = subprocess.run(f'{CLEAN_ENV} "cd {run} && sbatch --parsable --nice={a.nice} run.sbatch"', shell=True,
                           capture_output=True, text=True).stdout.strip()
        print(f"  -> job {j}")
        name = os.path.basename(run)
        lines = open(jobs_file).read().splitlines()
        lines = [f"{j} {name}" if ln.split()[1] == name else ln for ln in lines]
        open(jobs_file, "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
