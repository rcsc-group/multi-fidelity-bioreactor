"""Tests of the hydrodynamics-informed kLa law (zero compute): ln kLa(C* = 0.10) = c + s ln EDR.
(1) Does EDR carry information beyond (rpm, theta)? Correlation of the residuals of ln kLa and ln EDR after removing a quadratic in (rpm, theta) [L10] / (rpm) [L9].
(2) Cross-level: law calibrated at L9 (slope s9), applied at L10 with ONE level offset fitted on the L10 training runs: leave-one-out over the L10 theta-7 runs.
(3) Mixing time: ln dtmix_0.50 vs ln EDR at L8 (88 runs) and L9.
(4) FROZEN predictions for the theta-5.5 L10 kLa runs (not yet finished): law and rpm+theta regression, from the 13 L10 runs and the measured EDR of the theta-5.5 hydrodynamic runs.
Usage: uv run python scripts/physics_law_more.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.hydro_dataset_v2 import qois  # noqa: E402
from scripts.physics_law_kla_edr import table  # noqa: E402


def poly(rpm, th, deg2=True):
    cols = [np.ones_like(rpm), rpm, th]
    if deg2:
        cols += [rpm ** 2, rpm * th, th ** 2]
    return np.column_stack(cols)


def resid(y, A):
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    return y - A @ c


def main():
    T = table()
    t10, t9, t8 = T[10], T[9], T[8]
    # (1) partial correlation
    for L, a in ((10, t10), (9, t9)):
        rpm, th = a[:, 0], a[:, 1]
        A = poly(rpm, th, deg2=False) if L == 10 else np.column_stack([np.ones_like(rpm), rpm, rpm ** 2])
        r1, r2 = resid(np.log(a[:, 2]), A), resid(np.log(a[:, 3]), A)
        print(f"(1) L{L}: partial correlation of ln kLa and ln EDR given (rpm, theta) terms: {np.corrcoef(r1, r2)[0, 1]:+.2f} (n {len(a)})")
    # (2) cross-level
    m10 = t10[t10[:, 1] == 7.0]
    k9 = {round(r[0], 2): r for r in t9 if r[1] == 7.0}
    pts = [(r, k9[round(r[0], 2)]) for r in m10 if round(r[0], 2) in k9]
    X9 = np.log(t9[:, 3])
    s9 = np.polyfit(X9, np.log(t9[:, 2]), 1)[0]
    err = []
    for i, (r10, r9) in enumerate(pts):
        others = [p for j, p in enumerate(pts) if j != i]
        off = np.mean([np.log(a[2]) - np.log(b[2]) - s9 * (np.log(a[3]) - np.log(b[3])) for a, b in others])
        pred = np.log(r9[2]) + s9 * (np.log(r10[3]) - np.log(r9[3])) + off
        err.append(pred - np.log(r10[2]))
    base = []
    for i, (r10, r9) in enumerate(pts):
        off = np.mean([np.log(a[2]) - np.log(b[2]) for j, (a, b) in enumerate(pts) if j != i])
        base.append(np.log(r9[2]) + off - np.log(r10[2]))
    print(f"(2) L9 -> L10, theta 7, {len(pts)} points, slope from L9 {s9:.2f}: LOO RMS ln error law+offset {np.sqrt(np.mean(np.square(err))):.3f}; L9 x constant ratio {np.sqrt(np.mean(np.square(base))):.3f}")
    # (3) mixing time
    for L, name in ((8, "t6_l8_*"), (9, "fig9_l9_rpm*")):
        rows = []
        for d in sorted((ROOT / "runs").glob(name)) + (sorted((ROOT / "runs").glob("t5_l8_*")) + sorted((ROOT / "runs").glob("e8_l8_*")) if L == 8 else []):
            if "_ext" in d.name:
                continue
            f = d / "results.json"
            if not f.exists():
                continue
            v = json.loads(f.read_text()).get("dtmix_0.50")
            try:
                h = qois(d.name)
            except Exception:
                continue
            if v and v > 0 and h:
                rows.append((h["edr_mean_t"], v))
        a = np.array(rows)
        s, c0 = np.polyfit(np.log(a[:, 0]), np.log(a[:, 1]), 1)
        print(f"(3) L{L}: ln dtmix_0.50 = {c0:+.2f} {s:+.2f} ln EDR, r {np.corrcoef(np.log(a[:, 0]), np.log(a[:, 1]))[0, 1]:+.2f} (n {len(a)}); theory eps^(-1/3) gives slope -0.33")
    # (4) frozen predictions
    c = np.polyfit(np.log(t10[:, 3]), np.log(t10[:, 2]), 1)
    A = poly(t10[:, 0], t10[:, 1], deg2=False)
    cr, *_ = np.linalg.lstsq(A, np.log(t10[:, 2]), rcond=None)
    out = {}
    for r in (17.5, 20, 22.5, 27.5, 30, 32.5):
        run = f"x8_l10_rpm{r:g}_th5.5"
        if not (ROOT / "runs" / run / "results.json").exists():
            continue
        h = qois(run)
        out[f"{r:g}"] = dict(law=float(np.exp(np.polyval(c, np.log(h["edr_mean_t"])))), rpm_theta=float(np.exp(np.array([1, r, 5.5]) @ cr)), edr=float(h["edr_mean_t"]))
    (ROOT / "experiments/multifidelity/physics_law_frozen_theta55.json").write_text(json.dumps(dict(slope=float(c[0]), intercept=float(c[1]), n_train=len(t10), predictions=out), indent=1))
    print(f"(4) law from the {len(t10)} L10 runs: ln kLa = {c[1]:+.2f} + {c[0]:.2f} ln EDR; frozen predictions (1/h):", {k: (round(v['law'], 1), round(v['rpm_theta'], 1)) for k, v in out.items()})


if __name__ == "__main__":
    main()
