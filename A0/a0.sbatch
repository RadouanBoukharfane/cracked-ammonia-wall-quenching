#!/bin/bash
#SBATCH --job-name=A0
#SBATCH --account=msdaphyik-frgj4qkkjni-premium-cpu
#SBATCH --qos=premium-cpu
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=56
#SBATCH --exclusive
#SBATCH --time=12:00:00
#SBATCH --exclude=slurm-compute-h24b5-u30-svn3,slurm-compute-h21a5-u5-svn2,slurm-compute-h21c8-u5-svn4
#SBATCH --output=slurm-%x-%j.out
# usage: sbatch a0.sbatch [--no-soret] [--validate]
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=../.venv/bin/python
MECH=../mechanisms/mei2021/mei2021.yaml
SORET=""
for a in "$@"; do
  case $a in
    --validate) $PY run_a0.py validate --mech $MECH -j 56 ;;
    --no-soret) SORET="--no-soret" ;;
  esac
done
$PY run_a0.py run --mech $MECH --tag mei2021 $SORET -j 56
$PY run_a0.py retry --mech $MECH --tag mei2021 $SORET -j 56
