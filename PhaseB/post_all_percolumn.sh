#!/bin/bash
# Re-run post_hoq.py (per-column y_hrr) on all finished HOQ runs. Log: post_all_percolumn.log
cd "$(dirname "$0")"; export OMP_NUM_THREADS=1; PY=../.venv/bin/python
PLANAR=$(ls -d runs/HOQ_*planar* | grep -v "\.old")
UNST="runs/HOQ_g0.43_phi0.50_P10.0_T750_t20.0 runs/HOQ_g0.43_phi0.50_P10.0_T750_t26.0 runs/HOQ_g0.43_phi0.50_P10.0_T750_t32.0
runs/HOQ_g0.43_phi0.50_P10.0_T750_t38.0 runs/HOQ_g0.20_phi0.50_P10.0_T750_t20.0 runs/HOQ_g0.20_phi0.50_P10.0_T750_t26.0
runs/HOQ_g0.20_phi0.50_P10.0_T750_t38.0 runs/HOQ_g0.75_phi0.50_P10.0_T750_t26.0 runs/HOQ_g0.75_phi0.50_P10.0_T750_t38.0
runs/HOQ_g1.00_phi0.50_P10.0_T750_t26.0 runs/HOQ_g1.00_phi0.50_P10.0_T750_t38.0
runs/HOQ_g0.43_phi0.50_P10.0_T750_t26.0_Tw500 runs/HOQ_g0.43_phi0.50_P10.0_T750_t38.0_Tw500"
n=0; for r in $PLANAR; do $PY post_hoq.py $r 4 > $r/post_hoq.log 2>&1 & n=$((n+1)); [ $((n % 4)) -eq 0 ] && wait; done; wait
echo "planar done $(date)"
n=0; for r in $UNST; do $PY post_hoq.py $r 8 > $r/post_hoq.log 2>&1 & n=$((n+1)); [ $((n % 2)) -eq 0 ] && wait; done; wait
echo "unstable done $(date)"
