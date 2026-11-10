#!/bin/bash
#SBATCH -J pair-l8
#SBATCH -p batch
#SBATCH -n 1 -c 1
#SBATCH --mem=8G
#SBATCH -t 08:00:00
#SBATCH --exclude=node2356
#SBATCH -o /oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/logs/pair_l8_%j.out
# Driver for scripts/bo2d/pair_l8.py: advance every 10 min until all 12 L7 -> L8 pairs are done or Wed 2026-11-11 00:00.
cd /oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor
while [ $(date +%s) -lt $(date -d "2026-11-11 00:00" +%s) ]; do
  echo "== $(date -Is)"
  out=$(uv run python scripts/bo2d/pair_l8.py advance --go 2>&1 | grep -E "submitted|done|FAILED|no source|left|Error")
  echo "$out"
  echo "$out" | grep -q "^left 0$" && break
  sleep 600
done
echo "PAIR END $(date -Is)"
