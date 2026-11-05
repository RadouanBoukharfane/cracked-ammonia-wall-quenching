#!/bin/bash
# Waits for the 4 extra realisations (jobs_extra_realisations.txt), then post-processes them locally:
# post_hoq.py -> emissions_tq.py (adds the 4 runs to emissions_tq.csv) -> crest_at_quench.py. Log: post_extra_realisations.log
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
PY=../.venv/bin/python
RUNS=$(cut -d' ' -f1 jobs_extra_realisations.txt)
JOBS=$(cut -d' ' -f2 jobs_extra_realisations.txt | paste -sd,)
until [ -z "$(squeue -h -j $JOBS 2>/dev/null)" ]; do sleep 300; done
echo "jobs finished $(date)"
for r in $RUNS; do
  grep -q finalized $r/pele.log || echo "WARNING: $r did not reach the end (check before using)"
  $PY post_hoq.py $r 8 > $r/post_hoq.log 2>&1 &
done; wait
echo "post_hoq done $(date)"; for r in $RUNS; do head -1 $r/post_hoq.log; done
$PY emissions_tq.py 16 $RUNS 2>&1 | grep -v -i warn
$PY crest_at_quench.py 12 $RUNS 2>&1 | grep -v -i warn | grep "lead 2.0"
echo "all done $(date)"
