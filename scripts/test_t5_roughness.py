"""T5 roughness test: is dtmix_0.95(rpm) at L8, theta 7, real structure or effectively noise?

Pre-registration, diary 2026-10-08 (cont.), item T5, verbatim:
  "T5 (compute, ~100-300 core-h): roughness test at L8, same binary/protocol as kmix_l8. Dense rpm 20.625..31.875
  step 0.625 (15 new runs; with existing 20..32.5 at 2.5) + rpm perturbation 22.4/22.6/29.9/30.1 (4 runs). Reading
  fixed now:
    (a) REAL structure if the dense dtmix_0.95(rpm) curve has |second difference| at 0.625 spacing <= 3x the
    replicate log-sd (~0.02) except across <= 2 narrow features;
    (b) EFFECTIVELY NOISY if |log y(rpm +- 0.1) - log y(rpm)| or dense-neighbour jumps are >= 0.10 (comparable to
    the 2.5-rpm jumps 0.12-0.27)."

Data: theta 7, L8. Grid runs kmix_l8_rpm<r> (r = 15..37.5 step 2.5) and T5 runs t5_l8_rpm<r>_th7; one curve
(same protocol and binary). Dense curve = 20..32.5 at 0.625 (21 points: 6 grid + 15 T5).
Primary QoI dtmix_0.95. Reported only (not criteria): dtmix_0.50, dtmix_0.75, kLa_1T_10/25/50, and
tau_mean_t, tau_95_t, edr_mean_t via scripts.hydro_dataset_v2.qois (imported, unmodified).

Choices where the spec was silent:
  1. Replicate log-sd: dtmix_quality.json["0.95"]["replicate_sd"] entries L8@17.5, L8@25, L8@32.5; median used
     (0.0185); threshold = 3 x median.
  2. Second difference d2_i = log y_{i-1} - 2 log y_i + log y_{i+1} at the 19 interior dense points (no scaling).
  3. Exceedance = |d2| strictly > 3 x sd. A "narrow feature" = a maximal run of consecutive interior points that
     exceed (consecutive in the 0.625 index). Criterion (a) = number of such runs <= 2 (and exceedances all lie in them,
     which is then automatic). Zero exceedances also satisfies (a) (vacuous; flagged in the verdict text).
  4. Criterion (b): "noisy" = ANY of the pm0.1 pair differences (4 of them: 22.4-22.5, 22.6-22.5, 29.9-30,
     30.1-30) or dense-neighbour jumps |log y(r+0.625)-log y(r)| (20 of them) is >= 0.10; counts and offenders listed.
  5. NaN / non-positive values are dropped from that QoI's analysis (any pair or triple containing them is skipped).
  6. Hydro QoIs: qois() per run with its own params.json; tested for the 6 grid + 19 T5 runs; failures reported.
  7. Reported QoIs: same two metrics (second-difference exceedances vs 3 x the SAME dtmix_0.95 sd is not meaningful
     for other QoIs, so for them each is compared with 3 x its own L8 replicate sd only if available; replicates
     for these QoIs are not read, so only raw metric counts are given: n |d2| > 0.0555 (the dtmix_0.95 threshold),
     max |d2|, n pm0.1 and neighbour jumps >= 0.10, max jump).
  8. Figure: top row dtmix_0.95 with two zoom panels (22.4-22.6, 29.9-30.1); second row dtmix_0.50 and kLa_1T_10.
     Grid points larger markers, T5 smaller; pm0.1 points included in the curve line? No: the line joins the dense
     0.625 curve only; the pm0.1 points are drawn as separate markers.
Writes experiments/multifidelity/test_t5_roughness.{png,json}.
"""
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.figstyle import MK, GRID_KW, level_colour, rcparams, series_kw  # noqa: E402

GRID = [15 + 2.5 * i for i in range(10)]
DENSE = [20 + 0.625 * i for i in range(21)]
PERT = [22.4, 22.6, 29.9, 30.1]
QOIS = ["dtmix_0.95", "dtmix_0.50", "dtmix_0.75", "kLa_1T_10", "kLa_1T_25", "kLa_1T_50"]
HYD = ["tau_mean_t", "tau_95_t", "edr_mean_t"]
TH = 0.10


def run_id(r):
    return f"kmix_l8_rpm{r:g}" if any(abs(r - g) < 1e-9 for g in GRID) else f"t5_l8_rpm{r:g}_th7"


def is_grid(r):
    return any(abs(r - g) < 1e-9 for g in GRID)


def load(run):
    f = ROOT / "runs" / run / "results.json"
    return json.load(open(f)) if f.exists() else None


def metrics(vals, thr):
    """vals: dict rpm -> value (may be nan). Returns metric dict."""
    ly = {r: (math.log(v) if v == v and v > 0 else math.nan) for r, v in vals.items()}
    d2 = []
    for i in range(1, len(DENSE) - 1):
        a, b, c = (ly[DENSE[i + k]] for k in (-1, 0, 1))
        d2.append((DENSE[i], abs(a - 2 * b + c)))
    exc = [i for i, (_, d) in enumerate(d2) if d == d and d > thr]
    runs, cur = [], []
    for i in exc:
        if cur and i == cur[-1] + 1:
            cur.append(i)
        else:
            if cur:
                runs.append(cur)
            cur = [i]
    if cur:
        runs.append(cur)
    feats = [[d2[i][0] for i in rr] for rr in runs]
    pert = {}
    for p, ref in ((22.4, 22.5), (22.6, 22.5), (29.9, 30.0), (30.1, 30.0)):
        pert[f"{p:g}-{ref:g}"] = abs(ly[p] - ly[ref])
    jumps = {f"{DENSE[i]:g}->{DENSE[i+1]:g}": abs(ly[DENSE[i + 1]] - ly[DENSE[i]]) for i in range(len(DENSE) - 1)}
    big = {k: v for k, v in {**pert, **jumps}.items() if v == v and v >= TH}
    nanf = lambda dct: max([v for v in dct.values() if v == v], default=math.nan)
    return dict(d2={f"{r:g}": d for r, d in d2}, n_exceed=len(exc), exceed_rpm=[d2[i][0] for i in exc],
                features=feats, n_features=len(runs), max_d2=nanf(dict(d2)),
                pert=pert, jumps=jumps, n_b=len(big), b_offenders=big,
                max_pert=nanf(pert), max_jump=nanf(jumps),
                n_nan=sum(1 for v in ly.values() if v != v))


def main():
    sd_all = json.load(open(ROOT / "experiments/multifidelity/dtmix_quality.json"))["0.95"]["replicate_sd"]
    sds = {k: v for k, v in sd_all.items() if k.startswith("L8@")}
    sd = float(np.median(list(sds.values())))
    thr = 3 * sd
    print(f"L8 replicate log-sd (dtmix_0.95): {sds}; median {sd:.4f}; 3x = {thr:.4f}")
    rpms = DENSE + PERT
    res = {r: load(run_id(r)) for r in rpms}
    anomalies = [f"missing results.json: {run_id(r)}" for r in rpms if res[r] is None]
    data = {}
    for q in QOIS:
        data[q] = {r: (res[r].get(q, math.nan) if res[r] else math.nan) for r in rpms}
    # hydro
    hyd_err = []
    try:
        from scripts.hydro_dataset_v2 import qois
        for q in HYD:
            data[q] = {}
        for r in rpms:
            try:
                h = qois(run_id(r))
                for q in HYD:
                    data[q][r] = h[q]
            except Exception as e:
                hyd_err.append(f"{run_id(r)}: {type(e).__name__}: {e}")
                for q in HYD:
                    data[q][r] = math.nan
    except Exception as e:
        hyd_err.append(f"import failed: {e}")
        for q in HYD:
            data.pop(q, None)
    for q, d in data.items():
        for r, v in d.items():
            if not (v == v and v > 0):
                anomalies.append(f"{q} at {r:g} ({run_id(r)}): {v}")
    out = dict(sd_used=sd, sd_entries=sds, thr_a=thr, run_ids={f"{r:g}": run_id(r) for r in rpms},
               values={q: {f"{r:g}": v for r, v in d.items()} for q, d in data.items()},
               metrics={}, anomalies=anomalies, hydro_errors=hyd_err)
    for q in data:
        m = metrics(data[q], thr)
        out["metrics"][q] = m
        tag = "PRIMARY" if q == "dtmix_0.95" else "reported"
        print(f"\n[{q}] ({tag}) n|d2|>{thr:.4f}: {m['n_exceed']} of 19 at rpm {m['exceed_rpm']}; narrow features "
              f"{m['n_features']}: {m['features']}; max|d2| {m['max_d2']:.3f}")
        print(f"  pm0.1 pairs: " + ", ".join(f"{k}: {v:.3f}" for k, v in m["pert"].items()))
        print(f"  max dense-neighbour jump {m['max_jump']:.3f}; (b) offenders >= {TH}: {m['n_b']}: "
              f"{ {k: round(v, 3) for k, v in m['b_offenders'].items()} }; n_nan {m['n_nan']}")
    p = out["metrics"]["dtmix_0.95"]
    a = p["n_features"] <= 2
    b = p["n_b"] > 0
    out["criterion_a_real_structure"] = a
    out["criterion_b_effectively_noisy"] = b
    if a and not b:
        v = "REAL STRUCTURE (a holds, b does not)"
    elif b and not a:
        v = "EFFECTIVELY NOISY (b holds, a does not)"
    elif a and b:
        v = "BOTH (a) and (b) hold: pre-registered reading does not discriminate"
    else:
        v = "NEITHER (a) nor (b) holds"
    if a and p["n_exceed"] == 0:
        v += " [note: (a) holds vacuously, zero exceedances]"
    out["verdict"] = v
    print(f"\nanomalies: {anomalies}\nhydro errors: {hyd_err}")
    print(f"\nVERDICT: (a) = {a} ({p['n_features']} narrow features, {p['n_exceed']} exceedances); "
          f"(b) = {b} ({p['n_b']} values >= {TH}).\n  -> {v}")
    json.dump(out, open(ROOT / "experiments/multifidelity/test_t5_roughness.json", "w"), indent=1)
    figure(data, thr)


def figure(data, thr):
    plt.rcParams.update(rcparams())
    col = level_colour(8)
    fig = plt.figure(figsize=(11, 7.5))
    gs = fig.add_gridspec(2, 3, width_ratios=[2.2, 1, 1], hspace=0.55, wspace=0.7)

    def draw(ax, q, mk, ylab, xr=None):
        kw = series_kw(col, mk, stat="mean", ours=True)
        d = data[q]
        dr = [r for r in DENSE if (xr is None or xr[0] <= r <= xr[1])]
        ax.plot(dr, [d[r] for r in dr], color=col, ls=":", lw=kw["lw"])
        g = [r for r in dr if is_grid(r)]
        t = [r for r in dr if not is_grid(r)]
        ax.plot(g, [d[r] for r in g], ls="none", marker=mk, ms=9, mfc="w", mec=col, mew=1.3, label="grid (2.5 rpm)")
        ax.plot(t, [d[r] for r in t], ls="none", marker=mk, ms=5, mfc="w", mec=col, mew=1.1, label="T5 (0.625 rpm)")
        pp = [r for r in PERT if xr is None or xr[0] <= r <= xr[1]]
        ax.plot(pp, [d[r] for r in pp], ls="none", marker=mk, ms=6, mfc="w", mec="k", mew=1.1, label="rpm ± 0.1")
        ax.set_yscale("log")
        ax.set_xlabel("rpm (min$^{-1}$)")
        ax.set_ylabel(ylab)
        ax.grid(**GRID_KW)
    ax0 = fig.add_subplot(gs[0, 0])
    ylab = "$\\Delta t_{mix,0.95}$ (s)"
    draw(ax0, "dtmix_0.95", MK["chi_0.95"], ylab)
    ax0.set_title("$\\Delta t_{mix,0.95}$")
    ax0.legend(loc="upper left", bbox_to_anchor=(0, -0.2), ncol=3, frameon=False, fontsize=8)
    for j, xr in enumerate(((22.3, 22.7), (29.8, 30.2))):
        ax = fig.add_subplot(gs[0, 1 + j])
        draw(ax, "dtmix_0.95", MK["chi_0.95"], ylab, xr)
        ax.set_title(f"{xr[0]:g}–{xr[1]:g} rpm")
        ax.set_xticks([xr[0] + 0.1, xr[0] + 0.2, xr[0] + 0.3])
    ax1 = fig.add_subplot(gs[1, 0:2])
    draw(ax1, "dtmix_0.50", MK["chi_0.50"], "$\\Delta t_{mix,0.50}$ (s)")
    ax1.set_title("$\\Delta t_{mix,0.50}$")
    ax2 = fig.add_subplot(gs[1, 2])
    draw(ax2, "kLa_1T_10", MK["cstar_10"], "$k_La$ (h$^{-1}$)")
    ax2.set_title("$k_La$, $C^*=10\\%$")
    fig.savefig(ROOT / "experiments/multifidelity/test_t5_roughness.png", dpi=150, bbox_inches="tight")


if __name__ == "__main__":
    main()
