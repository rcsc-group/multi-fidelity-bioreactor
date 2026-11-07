"""Audit (2026-11-07): runs restarted with a theta change (theta_max_prev != theta_max) and whether they are rows of
the kLa / hydro databases used by the fits. The solver's phase sin(w_bio_st t) jumps at such a restart."""
import json, sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "campaign2d"))
import campaign as C  # noqa: E402
bad = []
for f in (C.ROOT / "runs").glob("*/params.json"):
    try:
        p = json.loads(f.read_text())
    except Exception:
        continue
    a = (p.get("theta_max_prev") or [0])[0]
    if a and abs(a - p["theta_max"][0]) > 1e-6:
        bad.append(f.parent.name)
kla = set(pd.read_csv(C.DB_KLA).run); hyd = {r["run"] for r in json.loads(C.DB_HYDRO.read_text())}
print("cross-theta restarts:", len(bad), "| in kLa db:", sorted(set(bad) & kla), "| in hydro db:", sorted(set(bad) & hyd))
