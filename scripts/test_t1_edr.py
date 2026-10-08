"""T1 (pre-registered, diary 2026-10-08): test-A protocol on the energy dissipation rate.

Pre-registration (diary.md, quoted):
  "T1 (zero compute): test A protocol unchanged, QoI = EDR (mean energy dissipation rate, same window
   as hydro_dataset_v2: last 3 cycles before release), LF L8/L9 -> HF L10 (9 rpm). PASS as test A: MF
   median rel RMSE below both HF-only GP and LF x ratio in >= 1 of the 2 LF cases at n_HF = 3 and not
   worse than 1.5x the best baseline in the other."
QoI edr_mean_t = time mean over the window of the domain-mean dissipation (shear_stress.dat column 9,
scale rho U^3/L, W/m^3), added to scripts/hydro_dataset_v2.py. Designs, methods, metric, coverage:
identical to scripts/test_mf_l10_holdout.py (run_design is imported, not copied); n_HF = 3 (7 designs)
and 4 (21 designs).

Choices where the spec was silent:
  * "the other case": the PASS rule is read as: at least one LF case is an outright MF win (below both
    baselines) AND the remaining case has MF median <= 1.5 x min(HF-only, LF x ratio) medians
    (trivially true if it is also a win).
  * EDR is modelled in raw units (no log transform), as test A did for tau.
  * Figure: median-MF-error n_HF=3 design per LF case (as test A).
Usage: uv run python scripts/test_t1_edr.py
"""
from __future__ import annotations

import itertools
import json
import sys
import time
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from scripts import figstyle as fs                                   # noqa: E402
from scripts.test_mf_l10_holdout import (load, run_design, pred, ENDS, LFS, METHODS, OUT)  # noqa: E402

Q = "edr_mean_t"


def main():
    t0 = time.time()
    warnings.filterwarnings("ignore")
    D = load()
    hrpm = np.array(sorted(D[10]))
    interior = [r for r in hrpm if r not in ENDS]
    yh_all = np.array([D[10][r][Q] for r in hrpm])
    res, models = {}, {}
    for lfl in LFS:
        xl = np.array(sorted(D[lfl]))
        yl = np.array([D[lfl][r][Q] for r in xl])
        for n in (3, 4):
            designs, fails = [], 0
            for mid in itertools.combinations(interior, n - 2):
                tr = list(ENDS) + list(mid)
                try:
                    o, mods = run_design(xl, yl, hrpm, yh_all, tr)
                except Exception as e:  # noqa: BLE001
                    fails += 1
                    print("  FAIL", lfl, tr, e)
                    continue
                designs.append(o)
                models[(lfl, n, tuple(tr))] = mods
            summ = {}
            for m in METHODS:
                v = np.array([d["rel"][m] for d in designs])
                summ[m] = dict(median=float(np.median(v)), min=float(v.min()), max=float(v.max()))
            for m in ("MF", "HF-only GP"):
                summ[m]["coverage95"] = sum(d["cov_hits"][m] for d in designs) / sum(d["n_test"] for d in designs)
            res[f"L{lfl}|n{n}"] = dict(summary=summ, designs=designs, n_fail=fails)

    wins, within = {}, {}
    for lfl in LFS:
        s = res[f"L{lfl}|n3"]["summary"]
        base = min(s["HF-only GP"]["median"], s["LF x ratio"]["median"])
        wins[f"L{lfl}"] = bool(s["MF"]["median"] < s["HF-only GP"]["median"] and s["MF"]["median"] < s["LF x ratio"]["median"])
        within[f"L{lfl}"] = bool(s["MF"]["median"] <= 1.5 * base)
    passed = bool(any(wins.values()) and all(within.values()))

    lines = [f"{'case':8s}{'nHF':>4s}  " + "  ".join(f"{m:>22s}" for m in METHODS) + "  cov MF/HF"]
    for key, r in res.items():
        lf, n = key.split("|")
        s = r["summary"]
        cells = "  ".join(f"{s[m]['median']:7.1%} [{s[m]['min']:5.1%},{s[m]['max']:6.1%}]" for m in METHODS)
        lines.append(f"{lf:8s}{n[1:]:>4s}  {cells}  {s['MF']['coverage95']:.0%}/{s['HF-only GP']['coverage95']:.0%}")
    table = "\n".join(lines)
    print(table)
    print("MF wins n=3:", wins, " within 1.5x best baseline:", within)
    print("PASS" if passed else "FAIL")

    plt.rcParams.update(fs.rcparams())
    f, axes = plt.subplots(1, 2, figsize=(8.0, 3.4), sharey=True)
    xg = np.linspace(15, 37.5, 200)
    rep = {}
    for j, lfl in enumerate(LFS):
        r = res[f"L{lfl}|n3"]
        order = np.argsort([d["rel"]["MF"] for d in r["designs"]])
        d = r["designs"][order[len(order) // 2]]
        tr = tuple(d["train"])
        rep[f"L{lfl}"] = dict(train=d["train"], mf_rel_rmse=d["rel"]["MF"])
        mf, sf, ratio = models[(lfl, 3, tr)]
        xl = np.array(sorted(D[lfl]))
        yl = np.array([D[lfl][x][Q] for x in xl])
        trm = np.isin(hrpm, tr)
        mu, sd = pred(mf, xg)
        mus, _ = pred(sf, xg)
        ax = axes[j]
        cl, ch, mk = fs.level_colour(lfl), fs.level_colour(10), fs.MK["ediss"]
        ax.fill_between(xg, mu - 1.96 * sd, mu + 1.96 * sd, color=ch, alpha=0.18, lw=0, label="MF 95%")
        ax.plot(xg, mu, color=ch, lw=1.4, ls="--", label="MF")
        ax.plot(xg, mus, color=fs.NEUTRAL, lw=1.2, ls="-.", label="HF-only GP")
        ax.plot(xg, ratio * np.interp(xg, xl, yl), color=fs.NEUTRAL, lw=1.2, ls=(0, (1, 3)), label="LF x ratio")
        ax.plot(xl, yl, **fs.series_kw(cl, mk, stat="mean", ours=True, ls="none"), label=f"L{lfl}")
        ax.plot(hrpm[trm], yh_all[trm], **fs.series_kw(ch, mk, stat="mean", ours=True, ls="none"), label="L10 training")
        ax.plot(hrpm[~trm], yh_all[~trm], **fs.series_kw(ch, mk, stat="mean", ours=True, ls="none", ms=3.0, alpha=0.6),
                label="L10 held out")
        ax.set_title(f"LF L{lfl}", fontsize=10)
        ax.set_xlabel("rocking speed (rpm)")
        ax.grid(**fs.GRID_KW)
        if j == 0:
            ax.set_ylabel(r"mean dissipation rate (W/m$^3$)")
        else:
            ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=8)
    f.tight_layout()
    f.savefig(OUT / "test_t1_edr.png", dpi=200, bbox_inches="tight")
    rt = time.time() - t0
    json.dump(dict(pass_=passed, mf_wins_n3=wins, within_1p5x_n3=within, representative_designs_n3=rep,
                   results=res, table=table, runtime_s=rt), open(OUT / "test_t1_edr.json", "w"), indent=1)
    print(f"runtime {rt:.1f} s")


if __name__ == "__main__":
    main()
