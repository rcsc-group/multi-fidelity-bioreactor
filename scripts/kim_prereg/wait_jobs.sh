#!/bin/bash
# Wait until none of the given Slurm job ids is in the queue. Usage: bash wait_jobs.sh ID [ID ...]
# Lists ALL of the user's jobs (squeue -j errors once a job is purged) and retries a failed squeue call,
# so a transient slurmctld error is never read as "queue empty".
while true; do
  out=$(squeue -h -u "$USER" -o %i 2>/dev/null); rc=$?
  if [ $rc -eq 0 ]; then
    left=0
    for id in "$@"; do grep -qx "$id" <<< "$out" && left=1; done
    [ $left -eq 0 ] && break
  fi
  sleep 10
done
echo "jobs $* left the queue"
