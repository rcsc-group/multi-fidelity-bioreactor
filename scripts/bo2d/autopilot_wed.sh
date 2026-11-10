#!/bin/bash
# BO autopilot (2026-11-08): every 30 min advance tasks; when an arm has no running/new task, propose its next round.
# Stops at the deadline (Wed 2026-11-11 00:00) or when STOP file exists. Runs as its own SLURM job (autopilot.sbatch); the nohup version died with the session 2026-11-08 16:31.
cd /oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor
LOG=logs/bo_autopilot.log
while [ ! -e scripts/bo2d/state/STOP ] && [ $(date +%s) -lt $(date -d "2026-11-11 00:00" +%s) ]; do
  echo "== $(date -Is)" >> $LOG
  uv run python scripts/bo2d/loop.py advance --go 2>&1 | grep -E "submitted|done|FAILED|no source" >> $LOG
  for arm in l9only chain; do
    n=$(uv run python scripts/bo2d/loop.py status 2>/dev/null | grep "^$arm " | grep -cE " (running|new) ")
    if [ "$n" -lt 4 ]; then   # asynchronous: fill free slots (2026-11-09)
      uv run python scripts/bo2d/loop.py propose $arm --go 2>&1 | grep -E "^$arm|propose|submitted|Error" >> $LOG
    fi
  done
  sleep 900
done
echo "AUTOPILOT END $(date -Is)" >> $LOG
