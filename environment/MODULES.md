# Cluster environment (Toubkal, UM6P)

- Solver build and runs: `module purge; module load foss/2024a` (GCC 13.3, OpenMPI 5), `ulimit -s unlimited`, 56 MPI ranks per node.
- Python 3.12 virtual environment with `requirements.txt` (Cantera, yt, contourpy, numpy, pandas, scipy, matplotlib).
- Paper: `module load texlive/20220321-GCC-11.2.0`. Movie: `module load FFmpeg/4.3.2-GCCcore-11.2.0`.
