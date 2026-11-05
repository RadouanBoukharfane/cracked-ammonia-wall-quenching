#!/bin/bash
#SBATCH --job-name=fq_analyses
#SBATCH --account=msdaphyik-frgj4qkkjni-premium-cpu
#SBATCH --qos=premium-cpu
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=56
#SBATCH --exclusive
#SBATCH --time=03:00:00
#SBATCH --exclude=slurm-compute-h24b5-u30-svn3,slurm-compute-h21a5-u5-svn2,slurm-compute-h21c8-u5-svn4,slurm-compute-h23a5-u13-svn2,slurm-compute-h23a5-u13-svn4,slurm-compute-h23a5-u14-svn1,slurm-compute-h23a5-u14-svn3,slurm-compute-h21a5-u14-svn3
#SBATCH --output=slurm-%x-%j.out
# TODO items 5, 6, 1, 4 (post-processing of the existing database, no new simulation). Outputs in analysis/.
# usage: sbatch run_analyses.sbatch      (rerun stability_database.py + matalon.py when wave 2 is complete)
cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=../.venv/bin/python
echo "== item 5: stability database $(date)"; $PY stability_database.py 2>&1 | grep -v -i warn
echo "== item 6: Matalon theory $(date)";     $PY matalon.py 2>&1 | grep -v -i warn | tail -3
echo "== item 1: local quenching events $(date)"; $PY local_events.py 28 2>&1 | grep -v -i warn
echo "== item 4: crest/trough emissions and budgets $(date)"; $PY conditioned_emissions.py 28 2>&1 | grep -v -i warn
echo "== done $(date)"
