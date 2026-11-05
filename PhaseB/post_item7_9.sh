#!/bin/bash
# Local post-processing of the TODO items 7 (12 planar HOQ q10/q50/q90) and 9 (2 planar HOQ, Zhang mechanism).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
PY=../.venv/bin/python
RUNS="$(ls -d runs/HOQ_crest_*_q[159]0_planar) $(ls -d runs/HOQ_mech_*)"
n=0
for r in $RUNS; do
  $PY post_hoq.py $r 4 > $r/post_hoq.log 2>&1 &
  n=$((n+1)); [ $((n % 4)) -eq 0 ] && wait
done; wait
echo "post_hoq done $(date)"; for r in $RUNS; do echo "$r: $(tail -1 $r/post_hoq.log)"; done
