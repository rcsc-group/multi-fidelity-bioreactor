#!/bin/bash
# Polls tau_chain_watchdog.py until all 13 L9+L10 tau-sweep conditions are
# done (self-healing broken L10 chain segments along the way), then exits 0.
# --mail-type=END on this job fires the "sweep finished" email; a crash in
# any individual segment already emails separately via its own --mail-type=FAIL.
#SBATCH --job-name=tau-sweep-sentinel
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=8-00:00:00
#SBATCH --output=logs/sentinel_%j.out
#SBATCH --error=logs/sentinel_%j.err

set -uo pipefail
cd /oscar/data/dharri15/eaguerov/Github/BioReactor3D/dev/rocking-bioreactor-2d

while true; do
    date
    uv run python scripts/tau_chain_watchdog.py
    status=$?
    if [ "$status" -eq 0 ]; then
        echo "All 13 L9/L10 tau-sweep conditions complete."
        exit 0
    fi
    sleep 1800
done
