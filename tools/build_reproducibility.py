"""Assemble the reproducibility package of the flame-quenching study from the working tree.

  python tools/build_reproducibility.py [OUT]        (default OUT = ./reproducibility)

The package mirrors the working-tree layout, so that every script runs unchanged from its root:
  code + inputs  (git-tracked)  : scripts, solver set-ups, mechanisms, case lists, per-run case.json/input.2d/run.sbatch,
                                  paper sources;
  derived data   (archive)      : per-run post-processed files (dispersion.csv, front.npz, hoq_*.npz/csv, ...), A0
                                  flame properties, analysis tables -> data_manifest.tsv (+ SHA-256) and data tarball.
Raw solver output (plotfiles, checkpoints) is not included. Absolute cluster paths are replaced by environment
variables (PELE_HOME, FQ_ROOT) in the copied code; run records (run.sbatch) keep their original form with the
project prefix replaced by <FQ_ROOT>.
"""
import glob
import hashlib
import os
import re
import shutil
import subprocess
import sys

SRC = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
OUT = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(SRC, "reproducibility"))
PROJ = "/srv/lustre01/project/msdaphyik-um6p-st-msda-frgj4qkkjni/users/radouan.boukharfane"
PELE = PROJ + "/solvers/PeleLMeX"
FQ = PROJ + "/flame-quenching"

CODE_GLOBS = [
    "A0/a0lib.py", "A0/run_a0.py", "A0/plot_a0.py", "A0/*.sbatch", "A0/JOURNAL.md",
    "PhaseA/*.py", "PhaseA/*.sbatch", "PhaseA/README.md", "PhaseA/template/*", "PhaseA/cases/*",
    "PhaseB0/*.py", "PhaseB0/*.sbatch", "PhaseB0/README.md", "PhaseB0/template/*", "PhaseB0/jobs*.txt",
    "PhaseB/*.py", "PhaseB/*.sh", "PhaseB/*.sbatch", "PhaseB/README.md", "PhaseB/jobs*.txt", "PhaseB/runs_*.txt",
    "PhaseB/build/*/GNUmakefile", "PhaseB/build/*/pelelmex_prob.*", "PhaseB/build/*/build.sbatch",
    "PhaseA/build/*/GNUmakefile", "PhaseA/build/*/build.sbatch",
    "PhaseB0/build/*/GNUmakefile", "PhaseB0/build/*/pelelmex_prob.*", "PhaseB0/build/*/build.sbatch",
    "analysis/*.py", "analysis/*.sbatch",
    "paper/main.tex", "paper/refs.bib", "paper/make_figs.py",
    "highlights/make_highlights.py", "highlights/make_movie.py",
    "mechanisms/reduce_for_pele.py", "mechanisms/mei2021/*", "mechanisms/mei2021/orig/*", "mechanisms/zhang2021/*",
    "watchdog.sh", "tools/build_reproducibility.py",
]
RUN_INPUTS = ["case.json", "input.2d", "run.sbatch"]
RUN_DATA = {
    "PhaseA/runs": ["pmf.dat", "dispersion.csv", "quality.json", "front.npz"],
    "PhaseB0/runs": ["pmf.dat", "b0_timeseries.csv", "local_mixture.csv", "local_mixture_stats.txt", "AC_History.dat"],
    "PhaseB/runs": ["pmf.dat", "hoq_timeseries.csv", "hoq_wallflux.npz", "hrr_total.csv", "crest_at_quench.csv",
                    "coldwall.npz", "post_hoq.log"],
}
DATA_GLOBS = ["A0/results/validation_mei2021.csv", "A0/results/mei2021/a0_summary.csv", "A0/results/mei2021/cases/*",
              "A0/results/mei2021/profiles/*", "A0/results/mei2021_nosoret/*.csv",
              "A0/results/mei2021_nosoret/cases/*", "PhaseB/*.csv", "analysis/*.csv",
              "PhaseB/crest/*/case.json", "PhaseB/crest/*/pmf.dat", "PhaseB/crest/*/input.2d", "PhaseB/crest/*/*.csv"]
TEXT_EXT = (".py", ".sh", ".sbatch", ".md", ".txt", ".tex", ".bib", ".json", ".2d", ".csv", "GNUmakefile", ".H", ".cpp",
            ".yaml", ".inp", ".dat")


def sanitize_code(text):
    """Replace cluster paths in copied code by environment-driven ones (valid Python / make / shell)."""
    # Python: MECH_ROOT / mechanism yaml definitions spread over adjacent string literals
    text = re.sub(r'\(\s*"' + re.escape(PROJ) + r'/solvers/"\s*\n\s*"PeleLMeX/Submodules/PelePhysics/Mechanisms"\s*\)',
                  'os.path.join(os.environ.get("PELE_HOME", "PeleLMeX"), "Submodules/PelePhysics/Mechanisms")', text)
    text = re.sub(r'\(\s*"' + re.escape(PROJ) + r'/solvers/PeleLMeX/"\s*\n\s*"Submodules/PelePhysics/Mechanisms/([^"]*)"\s*\)',
                  r'os.path.join(os.environ.get("PELE_HOME", "PeleLMeX"), "Submodules/PelePhysics/Mechanisms/\1")', text)
    text = re.sub(r'"' + re.escape(PELE) + r'/Submodules/PelePhysics/Mechanisms([^"]*)"',
                  r'os.path.join(os.environ.get("PELE_HOME", "PeleLMeX"), "Submodules/PelePhysics/Mechanisms\1")', text)
    text = re.sub(r'"' + re.escape(PROJ) + r'/solvers/PeleLMeX([^"]*)"',
                  r'(os.environ.get("PELE_HOME", "PeleLMeX") + "\1")', text)
    text = text.replace("PELE_HOME ?= " + PELE, "PELE_HOME ?= $(HOME)/PeleLMeX")
    text = text.replace(FQ, "${FQ_ROOT}").replace(PELE, "${PELE_HOME}")
    return text


def needs_os_import(text):
    return "os.environ" in text and not re.search(r"^import os\b|^import .*\bos\b", text, re.M)


def copy(rel, sanitize=False):
    src, dst = os.path.join(SRC, rel), os.path.join(OUT, rel)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if sanitize and rel.endswith(TEXT_EXT):
        t = open(src, errors="replace").read()
        t2 = sanitize_code(t) if not rel.endswith((".csv", ".dat", ".yaml", ".inp")) else t
        if rel.endswith(".py") and needs_os_import(t2):
            t2 = re.sub(r'^("""(?:.|\n)*?"""\n)', r"\1import os\n", t2, count=1)
        open(dst, "w").write(t2)
    else:
        shutil.copy2(src, dst)
    return dst


def main():
    if os.path.exists(OUT):
        sys.exit(f"{OUT} exists: remove it first")
    os.makedirs(OUT)
    code, data = [], []
    for g in CODE_GLOBS:
        for f in sorted(glob.glob(os.path.join(SRC, g))):
            if os.path.isfile(f):
                code.append(os.path.relpath(f, SRC))
    for base, files in RUN_DATA.items():
        for run in sorted(glob.glob(os.path.join(SRC, base, "*"))):
            if not os.path.isdir(run):
                continue
            for f in RUN_INPUTS:
                if os.path.exists(os.path.join(run, f)):
                    code.append(os.path.relpath(os.path.join(run, f), SRC))
            for f in files:
                if os.path.exists(os.path.join(run, f)):
                    data.append(os.path.relpath(os.path.join(run, f), SRC))
    for g in DATA_GLOBS:
        for f in sorted(glob.glob(os.path.join(SRC, g))):
            if os.path.isfile(f):
                data.append(os.path.relpath(f, SRC))
    code = sorted(set(code)); data = sorted(set(data) - set(code))
    for rel in code:
        copy(rel, sanitize=True)
    for rel in data:
        copy(rel)
    for d in ("analysis/figs", "paper/figs", "PhaseA/figs", "PhaseB0/figs", "PhaseB/figs", "highlights"):
        os.makedirs(os.path.join(OUT, d), exist_ok=True)
        open(os.path.join(OUT, d, ".gitkeep"), "w").close()
    # PelePhysics mechanisms generated with CEPTR (Mei2021NH3, Zhang2021NH3), and the local transport patch
    pp = os.path.join(PELE, "Submodules", "PelePhysics")
    for mech in ("Mei2021NH3", "Zhang2021NH3"):
        shutil.copytree(os.path.join(pp, "Mechanisms", mech), os.path.join(OUT, "solver", "PelePhysics_Mechanisms", mech))
    os.makedirs(os.path.join(OUT, "solver"), exist_ok=True)
    patch = subprocess.run(["git", "-C", pp, "diff"], capture_output=True, text=True).stdout
    open(os.path.join(OUT, "solver", "pelephysics_transport.patch"), "w").write(patch)
    ver = []
    for name, path in (("PeleLMeX", PELE), ("PelePhysics", pp), ("AMReX", os.path.join(pp, "Submodules", "amrex")),
                       ("SUNDIALS", os.path.join(pp, "Submodules", "sundials")),
                       ("AMReX-Hydro", os.path.join(PELE, "Submodules", "AMReX-Hydro"))):
        h = subprocess.run(["git", "-C", path, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
        url = subprocess.run(["git", "-C", path, "remote", "get-url", "origin"], capture_output=True, text=True).stdout.strip()
        ver.append(f"| {name} | `{h}` | {url} |")
    open(os.path.join(OUT, "solver", "VERSIONS.md"), "w").write(
        "# Solver versions used for all simulations\n\n| component | commit | repository |\n|---|---|---|\n"
        + "\n".join(ver) + "\n\nPelePhysics carries the local patch `pelephysics_transport.patch` (transport parameters)"
        " and the CEPTR-generated mechanisms in `PelePhysics_Mechanisms/` (copy them to "
        "`$PELE_HOME/Submodules/PelePhysics/Mechanisms/`).\n")
    # python environment
    os.makedirs(os.path.join(OUT, "environment"), exist_ok=True)
    freeze = subprocess.run([os.path.join(SRC, ".venv", "bin", "pip"), "freeze"], capture_output=True, text=True).stdout
    open(os.path.join(OUT, "environment", "requirements.txt"), "w").write(freeze)
    # documentation
    shutil.copy2(os.path.join(SRC, "tools", "repro_README.md"), os.path.join(OUT, "README.md"))
    open(os.path.join(OUT, "LICENSE"), "w").write(
        "MIT License (code)\n\nCopyright (c) 2026 Radouan Boukharfane\n\nPermission is hereby granted, free of charge, to any "
        "person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the "
        "Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, "
        "distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is "
        "furnished to do so, subject to the following conditions:\n\nThe above copyright notice and this permission "
        "notice shall be included in all copies or substantial portions of the Software.\n\nTHE SOFTWARE IS PROVIDED "
        "\"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF "
        "MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR "
        "COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT "
        "OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE "
        "SOFTWARE.\n\nData (post-processed results, inputs): Creative Commons Attribution 4.0 International (CC BY 4.0).\n")
    open(os.path.join(OUT, "CITATION.cff"), "w").write(
        "cff-version: 1.2.0\nmessage: Please cite the paper and this package.\ntitle: >-\n  Reproducibility package: pulsed "
        "wall heat loads, delayed ammonia burnout and cold-wall nitrous oxide in unstable partially cracked ammonia "
        "flames\nauthors:\n  - family-names: Boukharfane\n    given-names: Radouan\n    affiliation: Mohammed VI "
        "Polytechnic University (UM6P)\nversion: 1.0.0\ndate-released: 2026-11-05\nlicense: MIT\n")
    open(os.path.join(OUT, "environment", "MODULES.md"), "w").write(
        "# Cluster environment (Toubkal, UM6P)\n\n- Solver build and runs: `module purge; module load foss/2024a` "
        "(GCC 13.3, OpenMPI 5), `ulimit -s unlimited`, 56 MPI ranks per node.\n- Python 3.12 virtual environment with "
        "`requirements.txt` (Cantera, yt, contourpy, numpy, pandas, scipy, matplotlib).\n- Paper: `module load "
        "texlive/20220321-GCC-11.2.0`. Movie: `module load FFmpeg/4.3.2-GCCcore-11.2.0`.\n")
    # manifests
    with open(os.path.join(OUT, "data_manifest.tsv"), "w") as fh:
        fh.write("path\tbytes\tsha256\n")
        for rel in data:
            p = os.path.join(OUT, rel)
            fh.write(f"{rel}\t{os.path.getsize(p)}\t{hashlib.sha256(open(p, 'rb').read()).hexdigest()}\n")
    open(os.path.join(OUT, ".gitignore"), "w").write(
        "# derived data are distributed as an archive (see data_manifest.tsv), not in git\n"
        + "\n".join(sorted({"/" + r for r in data})) + "\n__pycache__/\n*.pyc\npaper/figs/\nhighlights/*.png\n"
        "highlights/*.pdf\nhighlights/*.mp4\n")
    remaining = subprocess.run(["grep", "-rl", "--include=*.py", PROJ, OUT], capture_output=True, text=True).stdout
    print(f"code/input files: {len(code)}, data files: {len(data)}")
    print("python files still holding cluster paths:", remaining.strip() or "none")


if __name__ == "__main__":
    main()
