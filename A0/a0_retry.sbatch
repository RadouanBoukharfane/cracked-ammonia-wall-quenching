#!/bin/bash
#SBATCH --job-name=A0_retry
#SBATCH --account=msdaphyik-frgj4qkkjni-premium-cpu
#SBATCH --qos=premium-cpu
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=56
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --exclude=slurm-compute-h24b5-u30-svn3,slurm-compute-h21a5-u5-svn2,slurm-compute-h21c8-u5-svn4
#SBATCH --output=slurm-%x-%j.out
# failed + missing cases of both grids, continuation in phi, 10 min max per case
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=../.venv/bin/python
MECH=../mechanisms/mei2021/mei2021.yaml
$PY run_a0.py retry --mech $MECH --tag mei2021 -j 28 --timeout 600 &
$PY run_a0.py retry --mech $MECH --tag mei2021 --no-soret -j 28 --timeout 600 &
wait
