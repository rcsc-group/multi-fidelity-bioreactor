#!/bin/bash
# BO autopilot (2026-11-08): every 30 min advance tasks; when an arm has no running/new task, propose its next round.
# Stops at the deadline (Tue 2026-11-10 18:00) or when STOP file exists. Runs inside interactive job 8361099.
cd /oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor
LOG=logs/bo_autopilot.log
while [ ! -e scripts/bo2d/state/STOP ] && [ $(date +%s) -lt $(date -d "2026-11-10 18:00" +%s) ]; do
  echo "== $(date -Is)" >> $LOG
  srun --jobid=8361099 --overlap -n1 -c1 --mem=6G uv run python scripts/bo2d/loop.py advance --go 2>&1 | grep -E "submitted|done|FAILED|no source" >> $LOG
  for arm in l9only chain; do
    n=$(uv run python scripts/bo2d/loop.py status 2>/dev/null | grep "^$arm " | grep -cE " (running|new) ")
    if [ "$n" -eq 0 ]; then
      srun --jobid=8361099 --overlap -n1 -c1 --mem=6G uv run python scripts/bo2d/loop.py propose $arm --go 2>&1 | grep -E "^$arm|propose|submitted|Error" >> $LOG
    fi
  done
  sleep 1800
done
echo "AUTOPILOT END $(date -Is)" >> $LOG
