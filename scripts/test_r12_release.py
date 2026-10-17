"""Test R12 (pre-registered, diary 2026-10-16): does releasing the tracer at cycle 12 instead of 80 change the L8 scalar QoIs?

d_i = log QoI(r12 run) - log QoI(T5 run at the same rpm), 12 rpm values, theta 7.
RELEASE-NEUTRAL: |mean d| <= 0.10 and sd(d) <= 1.5 sqrt(2) s_E0 for both primaries (kLa_1T_10 s_E0 0.095, dtmix_0.50 s_E0 0.104).
BIASED: |mean d| >= 0.15 for either primary. Else INCONCLUSIVE.
Usage: uv run python scripts/test_r12_release.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.hydro_dataset_v2 import qois  # noqa: E402
from scripts.submit_r12 import RPMS  # noqa: E402

PRIMARY = {"kLa_1T_10": 0.095, "dtmix_0.50": 0.104}
REPORTED = ("kLa_1T_25", "dtmix_0.95", "tau_95_t", "edr_mean_t")


def value(run, q):
    if q.endswith("_t"):
        h = qois(run)
        return h[q] if h else np.nan
    f = ROOT / "runs" / run / "results.json"
    return float(json.loads(f.read_text()).get(q, np.nan)) if f.exists() else np.nan


def main():
    ids = json.loads((ROOT / "experiments/multifidelity/test_t5_roughness.json").read_text())["run_ids"]
    out = {}
    verdicts = []
    for q in list(PRIMARY) + list(REPORTED):
        d = np.array([np.log(value(f"r12_l8_rpm{r:g}_th7", q)) - np.log(value(ids[f"{r:g}"], q)) for r in RPMS])
        ok = np.isfinite(d)
        m, s, n = float(d[ok].mean()), float(d[ok].std(ddof=1)), int(ok.sum())
        se = s / np.sqrt(n)
        out[q] = dict(mean=m, sd=s, n=n, se=se, d=d.tolist())
        tag = "PRIMARY" if q in PRIMARY else "reported"
        extra = ""
        if q in PRIMARY:
            bar_sd = 1.5 * np.sqrt(2) * PRIMARY[q]
            neutral = abs(m) <= 0.10 and s <= bar_sd
            biased = abs(m) >= 0.15
            verdicts.append("BIASED" if biased else "NEUTRAL" if neutral else "INCONCLUSIVE")
            extra = f"  sd bar {bar_sd:.3f} -> {verdicts[-1]}"
        print(f"{q:11s} [{tag}] n {n}  mean d {m:+.3f} (se {se:.3f})  sd d {s:.3f}{extra}")
    final = "BIASED" if "BIASED" in verdicts else "RELEASE-NEUTRAL" if all(v == "NEUTRAL" for v in verdicts) else "INCONCLUSIVE"
    out["verdict"] = final
    print("R12 verdict:", final)
    (ROOT / "experiments/multifidelity/test_r12_release.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
