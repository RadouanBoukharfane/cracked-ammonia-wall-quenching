#!/bin/bash
#SBATCH --job-name=HOQ_post
#SBATCH --account=msdaphyik-frgj4qkkjni-premium-cpu
#SBATCH --qos=premium-cpu
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=56
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --exclude=slurm-compute-h24b5-u30-svn3,slurm-compute-h21a5-u5-svn2,slurm-compute-h21c8-u5-svn4,slurm-compute-h23a5-u13-svn2,slurm-compute-h23a5-u13-svn4,slurm-compute-h23a5-u14-svn1,slurm-compute-h23a5-u14-svn3,slurm-compute-h21a5-u14-svn3
#SBATCH --output=slurm-%x-%j.out
# usage: sbatch post_hoq.sbatch RUN_DIR [RUN_DIR ...]   (one after the other, 28 processes each)
cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=1
for r in "$@"; do
  ../.venv/bin/python post_hoq.py "$r" 28 > "$r/post_hoq.log" 2>&1; echo "$r exit=$?"; head -1 "$r/post_hoq.log"
done
