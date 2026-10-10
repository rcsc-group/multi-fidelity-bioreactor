"""Is the mean energy dissipation rate (EDR, W/m^3 = power per volume) a usable proxy for kLa across rpm?

Pre-registered (diary 2026-10-10): proxy usable if Spearman(EDR, kLa_1T_c) >= 0.8 at BOTH L8 and L9 for
c = 10 and 25 (same 10 rpm, theta 7, kmix_l8 / fig9_l9 runs). Also reported: log-log slope of kLa on EDR.
EDR comes from experiments/multifidelity/hydro_dataset_v2.json (last 3 cycles before release); kLa from
runs/<run>/results.json.

Usage: uv run python scripts/diag_edr_kla_proxy.py
"""
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
SRC = {8: "kmix_l8", 9: "fig9_l9"}
KLA = ("kLa_1T_10", "kLa_1T_25", "kLa_1T_50")


def main():
    hd = json.load(open(ROOT / "experiments/multifidelity/hydro_dataset_v2.json"))
    out = {}
    for lev, prefix in SRC.items():
        rows = sorted((r for r in hd if r["level"] == lev and r["run"].startswith(prefix)), key=lambda r: r["rpm"])
        edr = np.array([r["edr_mean_t"] for r in rows])
        for k in KLA:
            y = np.array([json.load(open(ROOT / "runs" / r["run"] / "results.json")).get(k, np.nan) for r in rows],
                         float)
            m = np.isfinite(y)
            rho = spearmanr(edr[m], y[m]).statistic
            slope = np.polyfit(np.log(edr[m]), np.log(y[m]), 1)[0]
            out[f"L{lev}|{k}"] = dict(n=int(m.sum()), spearman=float(rho), loglog_slope=float(slope))
            print(f"L{lev} {k:10s} n={m.sum():2d} Spearman {rho:+.2f}  log-log slope {slope:+.2f}")
    ok = all(out[f"L{lev}|{k}"]["spearman"] >= 0.8 for lev in SRC for k in KLA[:2])
    print("PASS" if ok else "FAIL")
    out["pass"] = ok
    json.dump(out, open(ROOT / "experiments/multifidelity/diag_edr_kla_proxy.json", "w"), indent=1)


if __name__ == "__main__":
    main()
