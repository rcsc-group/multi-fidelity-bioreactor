"""Join diagnostics for a continuation: compare run B (joined with its parent) with the uninterrupted reference.

Prints QoIs of both, relative differences, and C*, chi, sigma^2 around the join, plus the first rows of B's
tr_oxy.dat vs the last of A's, and the reference's values at the same times.
Usage: uv run python scripts/diag_continue_join.py REF B
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
from scripts import postprocess as pp  # noqa: E402

KEYS = ["dtmix_0.50", "dtmix_0.75", "dtmix_0.95", "kLa_1T_10", "kLa_1T_25", "kLa_1T_50", "kLa_10", "kLa_25",
        "kLa_50", "tau_95_qss", "tau_95_max", "vor_streaming", "ediss_mean_qss", "vel_rms_qss"]


def series(d: Path):
    t, c2, c2s = pp._joined_cols(d, "tr_oxy.dat", [pp._COL_C2_LIQ_SUM, pp._COL_C2_LIQ_SUM2])
    _, ox = pp._joined_cols(d, "tr_oxy.dat", [pp._COL_OXY_LIQ])
    _, f = pp._joined_cols(d, "vol_frac_interf.dat", [pp._COL_F_LIQ])
    fm = f.mean()
    s2 = c2s / fm - (c2 / fm) ** 2
    i0 = np.where(c2 > 1e-10 * fm)[0][0]
    return t, ox / fm, 1 - s2 / s2[i0], s2


def main() -> None:
    ref, b = (ROOT / "runs" / a for a in sys.argv[1:3])
    rr, rb = (json.loads((d / "results.json").read_text()) for d in (ref, b))
    print(f"{'qoi':16s}{'ref':>14s}{'B':>14s}{'rel diff %':>12s}")
    for k in KEYS:
        x, y = rr.get(k), rb.get(k)
        print(f"{k:16s}{x:14.6g}{y:14.6g}{100 * (y - x) / x:12.2f}" if x else f"{k:16s}{x} {y}")
    tr, Cr, chir, s2r = series(ref)
    tb, Cb, chib, s2b = series(b)
    seg = pp._segment_chain(b)
    tj = float(pp._load_col(seg[-1] / "tr_oxy.dat", pp._COL_T)[0])
    print(f"\njoin at t = {tj:.6f}; segments {[s.name for s in seg]}")
    print(f"{'t':>10s}{'C*_ref':>11s}{'C*_B':>11s}{'chi_ref':>11s}{'chi_B':>11s}{'s2_ref':>12s}{'s2_B':>12s}")
    ib = np.searchsorted(tb, tj)
    for k in range(ib - 4, ib + 5):
        j = int(np.argmin(abs(tr - tb[k])))
        print(f"{tb[k]:10.4f}{Cr[j]:11.6f}{Cb[k]:11.6f}{chir[j]:11.6f}{chib[k]:11.6f}{s2r[j]:12.6e}{s2b[k]:12.6e}")
    # growth of the C*, chi differences after the join
    print("\nafter join: max |dC*|, max |dchi| over successive windows (B - ref on the same t)")
    jj = np.array([int(np.argmin(abs(tr - x))) for x in tb])
    for lo, hi in [(tj, tj + 2), (tj + 2, tj + 10), (tj + 10, tj + 40), (tj + 40, 1e9)]:
        m = (tb >= lo) & (tb < hi)
        if m.any():
            print(f"  t in [{lo:.1f},{min(hi, tb[-1]):.1f}): dC* mean {np.mean(Cb[m] - Cr[jj[m]]):+.5f} "
                  f"max {np.max(abs(Cb[m] - Cr[jj[m]])):.5f}; dchi mean {np.mean(chib[m] - chir[jj[m]]):+.5f} "
                  f"max {np.max(abs(chib[m] - chir[jj[m]])):.5f}")
    print("\nfirst rows of B tr_oxy.dat:")
    for line in (seg[-1] / "tr_oxy.dat").read_text().splitlines()[:4]:
        print("  ", line[:110])
    print("last rows of A tr_oxy.dat:")
    for line in (seg[-2] / "tr_oxy.dat").read_text().splitlines()[-2:]:
        print("  ", line[:110])


if __name__ == "__main__":
    main()
