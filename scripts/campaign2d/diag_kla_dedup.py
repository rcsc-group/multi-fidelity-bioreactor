"""Diagnostic (2026-11-04): why kLa_1T_10 max 2 sigma/mean moved from 0.76 to 2.40 when the six duplicate rows were
dropped. Fits the same model with and without dedup; prints the data layout, the nearest data to the corner
(rpm 15, theta 2), the model hyperparameters where exposed, and crit at a few points. Compute node only.
Usage: uv run python scripts/campaign2d/diag_kla_dedup.py
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import campaign as c  # noqa: E402

QOI = "kLa_1T_10"
PTS = np.array([[15.0, 2.0], [15.0, 4.0], [22.5, 2.0], [37.5, 7.0], [30.0, 5.0]])


def hyper(m):
    out = {}
    for name in dir(m):
        low = name.lower()
        if name.startswith("_") or not any(k in low for k in ("len", "scale", "noise", "sig", "kern", "hyp", "amp")):
            continue
        try:
            v = getattr(m, name)
        except Exception:
            continue
        if isinstance(v, (int, float, np.floating)) or (isinstance(v, np.ndarray) and v.size <= 8):
            out[name] = np.round(np.asarray(v, float), 4).tolist()
    return out


def run(label, d):
    print(f"\n== {label}: rows {len(d)}, runs {d.run.nunique()}")
    print("   by level:", d.groupby("level").size().to_dict(), " theta<4 rows:", int((d.theta < 4).sum()))
    print("   L8 theta values:", sorted(d[d.level == 8].theta.round(2).unique().tolist()))
    dist = np.hypot((d.rpm - 15.0) / 22.5, (d.theta - 2.0) / 5.0)
    k = int(np.argmin(dist.values))
    print(f"   nearest data to (15, 2): rpm {d.rpm.iloc[k]}, theta {d.theta.iloc[k]}, level {d.level.iloc[k]}, "
          f"normalised distance {dist.iloc[k]:.3f}")
    m = c.fit_qoi(d, QOI)
    lf = next((getattr(m, a) for a in dir(m) if a.lower().startswith("lf")
               and hasattr(getattr(m, a), "kernel")), None)
    if lf is not None:
        kv = {}
        for k, v in vars(lf.kernel).items():
            try:
                kv[k] = np.round(np.asarray(v, float), 4).tolist()
            except (TypeError, ValueError):
                pass
        print("   LF kernel params:", kv)
        print("   LF hyperparams (noise, y_std):",
              " LF noise_std:", np.round(float(np.asarray(lf.noise_std).ravel()[0]), 5),
              " LF y_std:", np.round(float(np.asarray(getattr(lf, "y_std", np.nan)).ravel()[0]), 5))
    else:
        print("   LF model attribute not found; attributes:", [a for a in dir(m) if not a.startswith("_")][:40])
    cq = c.crit(m, PTS)
    for p, v in zip(PTS, cq):
        print(f"   crit at rpm {p[0]:5.1f} theta {p[1]:4.1f}: {v:.3f}")
    r = c.max_crit(m)
    print(f"   continuous max {r['max']:.3f} at rpm {r['rpm']:.2f} theta {r['theta']:.2f}")


def main():
    import pandas as pd
    kla = pd.read_csv(c.DB_KLA)
    raw = kla[["level", "rpm", "theta", "run", QOI]].rename(columns={QOI: "y"}).dropna()
    raw = raw[raw.level.isin([8, 9, 10]) & (raw.rpm >= c.RPM[0]) & (raw.rpm <= c.RPM[1])
              & (raw.theta >= c.THETA[0]) & (raw.theta <= c.THETA[1]) & (raw.y > 0)].reset_index(drop=True)
    run("WITH duplicates (before fix)", raw)
    run("DEDUP (current)", c.dedup_runs(raw))


if __name__ == "__main__":
    main()
