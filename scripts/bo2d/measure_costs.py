"""Measured core-h per level from sacct (2026-11-08), replacing the guessed costs (L9 guessed 77, measured ~4x more)."""
import json
import subprocess
from pathlib import Path
import numpy as np
jobs = [json.loads(l) for l in Path("scripts/bo2d/state/jobs.jsonl").read_text().splitlines() if l.strip()]
out = {}
for j in jobs:
    r = subprocess.run(["sacct", "-X", "-n", "-P", "-j", str(j.get("job")), "-o", "State,ElapsedRaw,NCPUS"],
                       capture_output=True, text=True).stdout.split()
    if not r or not r[0].startswith("COMPLETED"):
        continue
    st, el, nc = r[0].split("|")
    key = (j.get("level"), "spin" if j.get("spin_only") else ("warm" if j.get("source") else "cold"))
    out.setdefault(key, []).append(int(el) * int(nc) / 3600)
for k, v in sorted(out.items(), key=str):
    print(k, f"n={len(v)} median {np.median(v):.0f} core-h  range {min(v):.0f}-{max(v):.0f}")
