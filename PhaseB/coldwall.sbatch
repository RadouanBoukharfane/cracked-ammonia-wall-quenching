#!/bin/bash
#SBATCH --job-name=coldwall
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
# usage: sbatch coldwall.sbatch   (post_hoq for the extended run, then coldwall_analysis on every planar T_u = 750 K run)
cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=1

../.venv/bin/python coldwall_analysis.py 56 $(ls -d runs/HOQ_g*_phi0.50_P10.0_T750_planar*) 2>&1 | grep -v -i warn
