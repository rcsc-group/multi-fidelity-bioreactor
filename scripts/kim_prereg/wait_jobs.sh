#!/bin/bash
# Wait until none of the given Slurm job ids is in the queue. Usage: bash wait_jobs.sh ID [ID ...]
ids=$(IFS=,; echo "$*")
while [ -n "$(squeue -h -j "$ids" -o %i 2>/dev/null)" ]; do sleep 10; done
echo "jobs $ids left the queue"
