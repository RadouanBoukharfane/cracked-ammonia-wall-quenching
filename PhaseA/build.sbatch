#!/bin/bash
#SBATCH --job-name=A_build
#SBATCH --account=msdaphyik-frgj4qkkjni-premium-cpu
#SBATCH --qos=premium-cpu
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=56
#SBATCH --exclusive
#SBATCH --time=03:00:00
#SBATCH --exclude=slurm-compute-h24b5-u30-svn3,slurm-compute-h21a5-u5-svn2,slurm-compute-h21c8-u5-svn4
#SBATCH --output=slurm-%x-%j.out
# Builds one PeleLMeX executable per pressure (PLOG reactions) in build/<PDIR>/ -> bin/
module purge; module load foss/2024a
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"; mkdir -p bin; date
# third-party libs (SUNDIALS) once, in the first build dir
P0=1_000atm
mkdir -p build/$P0 && cp template/{GNUmakefile,pelelmex_prob.H,pelelmex_prob.cpp} build/$P0/
(cd build/$P0 && make TPL COMP=gnu PDIR=$P0 > build_tpl.log 2>&1) || { echo TPL failed; exit 1; }
for P in 1_000atm 5_000atm 10_000atm 0_987atm 9_869atm; do
  mkdir -p build/$P && cp template/{GNUmakefile,pelelmex_prob.H,pelelmex_prob.cpp} build/$P/
  (cd build/$P && make COMP=gnu PDIR=$P -j 11 > build.log 2>&1 && cp PeleLMeX2d.gnu.MPI.Mei2021NH3_$P.ex ../../bin/ && echo "$P ok" || echo "$P FAILED") &
done
wait; ls -la bin; date
