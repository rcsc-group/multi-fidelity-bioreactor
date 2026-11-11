"""Is there an L9 kLa landscape off theta = 7? (2026-11-10 23:00). All L9 runs with finite kLa_1T_10 (bo1_l9, bo4_l9,
fig9_l9 theta = 7 line). Per set: sd ln kLa vs the noise floor (campaign.noise_floor); Spearman of kLa with rpm, theta
and with the SAME run's hydro features. Same for dtmix_0.25 (the rate 1/t) for contrast.
  uv run python scripts/bo2d/diag_l9_kla_structure.py
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "scripts" / "campaign2d")]
import campaign as C  # noqa: E402

FE = ["vel_rms_qss", "edr_95", "ediss_mean_qss", "tau_95_qss"]
rows = []
for g in ("bo1_l9_*", "bo4_l9_*", "fig9_l9_rpm*"):
    for d in sorted((ROOT / "runs").glob(g)):
        if "failed" in d.name or not (d / "results.json").exists():
            continue
        r = json.loads((d / "results.json").read_text()); p = json.loads((d / "params.json").read_text())
        rpm = p["omega_b"] * 60 / (2 * np.pi); th = p["theta_max"][0]
        mix = C.mix_value(d, "dtmix_0.25")
        rows.append(dict(set=g.split("_")[0], rpm=rpm, th=th, kla=r.get("kLa_1T_10"), mix=mix, **{k: r.get(k) for k in FE}))
for q, fl in (("kla", C.noise_floor("kLa_1T_10")), ("mix", C.noise_floor("dtmix_0.25"))):
    for sel, f in (("off-line", lambda r: abs(r["th"] - 7) > 0.01), ("theta=7", lambda r: abs(r["th"] - 7) <= 0.01)):
        R = [r for r in rows if f(r) and r[q] is not None and np.isfinite(r[q]) and r[q] > 0]
        if len(R) < 5:
            continue
        y = np.log([r[q] for r in R])
        s = f"{q} {sel:8s} n={len(R):2d} sd ln {np.std(y):.2f} (floor {fl:.3f})"
        for k in ["rpm", "th"] + FE:
            v = np.array([r[k] if r[k] is not None else np.nan for r in R], float); ok = np.isfinite(v)
            s += f" | {k} {spearmanr(v[ok], y[ok]).correlation:+.2f}"
        print(s)
