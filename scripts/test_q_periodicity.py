"""Test Q (periodicity per rpm), pre-registered in diary.md.

Metric: e(t) = ||U(t+T) - U(t)|| / ||U(t)||, U = (ux, uy) over cells, U(t+T)
by linear interpolation in time between the two snapshots bracketing t+T.
Liquid version: only cells with f > 0.5 and cs > 0.5 at t AND at both snapshots
bracketing t+T.  All-cell version also reported.  Median over every phase t
whose t+T lies inside the record.
Floor: ||interp(t-dt, t+dt) - U(t)|| / ||U(t)|| (liquid in all three), median
over snapshots.  If floor > 1 %, PERIODIC bar = 2 x floor (printed).
Class: PERIODIC if median e < bar; APERIODIC if > 10 %; else QUASI-PERIODIC.

T is the period in the snapshot (nondimensional) time: pp._t_scales(params)
-> T_per_nd (omega_b alone gives the DIMENSIONAL period in s).

Usage:
  uv run python scripts/test_q_periodicity.py                    # 10 flow_l8 runs
  uv run python scripts/test_q_periodicity.py --runs rec_l8_rpm22.4_th7 --periods 4
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

OUT = ROOT / "experiments/multifidelity"
BAR_PERIODIC, BAR_APERIODIC = 0.01, 0.10


# ------------------------------------------------------------- pure metric
def _norm(x):
    return np.sqrt(np.sum(x * x))


def periodicity_error(times, U, mask, T):
    """times (nt,), U (nt, 2, ny, nx), mask (nt, ny, nx) bool, period T.
    Returns (phases t, e_liquid, e_all) for every t with t+T <= times[-1]."""
    times = np.asarray(times, float)
    ph, el, ea = [], [], []
    for i, t in enumerate(times):
        tt = t + T
        if tt > times[-1] + 1e-12:
            break
        j = int(np.searchsorted(times, tt - 1e-12))
        j = min(max(j, 1), len(times) - 1)
        w = (tt - times[j - 1]) / (times[j] - times[j - 1])
        w = min(max(w, 0.0), 1.0)
        Up = (1 - w) * U[j - 1] + w * U[j]
        d = Up - U[i]
        m = mask[i] & mask[j - 1] & mask[j]
        ph.append(t)
        ea.append(_norm(d) / _norm(U[i]))
        den = _norm(U[i][:, m])
        el.append(_norm(d[:, m]) / den if m.any() and den > 0 else np.nan)
    return np.array(ph), np.array(el), np.array(ea)


def floor_error(times, U, mask):
    """Error of rebuilding snapshot k from k-1 and k+1 (liquid cells in all
    three), for k = 1..nt-2."""
    times = np.asarray(times, float)
    out = []
    for k in range(1, len(times) - 1):
        w = (times[k] - times[k - 1]) / (times[k + 1] - times[k - 1])
        Ui = (1 - w) * U[k - 1] + w * U[k + 1]
        m = mask[k - 1] & mask[k] & mask[k + 1]
        den = _norm(U[k][:, m])
        out.append(_norm((Ui - U[k])[:, m]) / den if m.any() and den > 0 else np.nan)
    return np.array(out)


def bar_from_floor(floor):
    """(PERIODIC bar, raised?)"""
    return (2 * floor, True) if floor > BAR_PERIODIC else (BAR_PERIODIC, False)


def classify(e_med, floor):
    bar, _ = bar_from_floor(floor)
    if e_med < bar:
        return "PERIODIC", bar
    if e_med > BAR_APERIODIC:
        return "APERIODIC", bar
    return "QUASI-PERIODIC", bar


# -------------------------------------------------------------- run access
def load_run(run_id, periods=None):
    from scripts import chi_scale as cs_
    from scripts import postprocess as pp
    loc = cs_.locate(run_id)
    if loc["fields_dir"] is None or loc["data_dir"] is None:
        return None
    params = json.loads((loc["data_dir"] / "params.json").read_text())
    _, T = pp._t_scales(params)
    files = sorted(Path(loc["fields_dir"]).glob("snap_*.bin"))
    U, M, ts = [], [], []
    t0 = None
    for p in files:
        s = cs_.read_snapshot(p, ["ux", "uy", "f", "cs"])
        if t0 is None:
            t0 = s["t"]
        if periods is not None and s["t"] > t0 + periods * T + 1e-9:
            break
        ts.append(s["t"])
        U.append(np.stack([s["ux"], s["uy"]]).astype(np.float32))
        M.append((s["f"] > 0.5) & (s["cs"] > 0.5))
    ts = np.array(ts)
    return {"params": params, "T": float(T), "t": ts, "U": np.array(U), "mask": np.array(M),
            "rpm": params["omega_b"] * 60 / (2 * np.pi),
            "fpp": float(T / np.median(np.diff(ts)))}


def analyse(d):
    ph, e, e_all = periodicity_error(d["t"], d["U"], d["mask"], d["T"])
    fl = floor_error(d["t"], d["U"], d["mask"])
    fl_all = floor_error(d["t"], d["U"], np.ones_like(d["mask"]))
    em, floor = float(np.nanmedian(e)), float(np.nanmedian(fl))
    cls, bar = classify(em, floor)
    q1, q3 = np.nanpercentile(e, [25, 75])
    return {"rpm": d["rpm"], "T_nd": d["T"], "frames_per_period": d["fpp"], "n_snap": len(d["t"]),
            "median_e": em, "q25_e": float(q1), "q75_e": float(q3),
            "median_e_allcells": float(np.nanmedian(e_all)),
            "floor": floor, "floor_allcells": float(np.nanmedian(fl_all)),
            "bar": bar, "bar_raised": bool(bar > BAR_PERIODIC), "class": cls,
            "phase_t": ph.tolist(), "e_liquid": e.tolist(), "e_all": e_all.tolist()}


def make_figure(res, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scripts import figstyle as fs
    plt.rcParams.update(fs.rcparams())
    rows = sorted(res.values(), key=lambda r: r["rpm"])
    x = np.array([r["rpm"] for r in rows])
    med = np.array([r["median_e"] for r in rows])
    lo = med - np.array([r["q25_e"] for r in rows])
    hi = np.array([r["q75_e"] for r in rows]) - med
    c = fs.level_colour(8)
    f, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.errorbar(x, med, yerr=[lo, hi], capsize=3, label="median $e$ (IQR)",
                **fs.series_kw(c, fs.MK["l2"], stat="mean", ours=True))
    ax.plot(x, [r["floor"] for r in rows], label="interpolation floor",
            **fs.series_kw(c, "+", stat="max", ours=True, mfc=c))
    ax.axhline(BAR_PERIODIC, color=fs.NEUTRAL, lw=0.9, ls="--", label="1%")
    ax.axhline(BAR_APERIODIC, color=fs.NEUTRAL, lw=0.9, ls="-.", label="10%")
    ax.set_yscale("log")
    ax.set_xlabel("rocking rate (rpm)")
    ax.set_ylabel(r"$\|U(t+T)-U(t)\|\,/\,\|U(t)\|$")
    ax.grid(**fs.GRID_KW)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    path.parent.mkdir(parents=True, exist_ok=True)
    f.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(f)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="*", default=None)
    ap.add_argument("--periods", type=float, default=None,
                    help="use only the first N periods of the record (rec_* runs hold 12)")
    ap.add_argument("--tag", default="test_q_periodicity")
    a = ap.parse_args()
    runs = a.runs or sorted(p.name for p in (ROOT / "runs").glob("flow_l8_rpm*_th7"))
    res, missing = {}, []
    for r in runs:
        d = load_run(r, a.periods)
        if d is None:
            missing.append(r)
            continue
        res[r] = analyse(d)
    print(f"{'run':28s} {'rpm':>6s} {'fpp':>6s} {'med e':>9s} {'IQR e':>19s} {'e_allcell':>9s} "
          f"{'floor':>9s} {'bar':>8s}  class")
    for r, v in res.items():
        print(f"{r:28s} {v['rpm']:6.2f} {v['frames_per_period']:6.1f} {v['median_e']:9.3e} "
              f"[{v['q25_e']:8.2e},{v['q75_e']:8.2e}] {v['median_e_allcells']:9.3e} "
              f"{v['floor']:9.3e} {v['bar']:8.2e}  {v['class']}")
        if v["bar_raised"]:
            print(f"  NOTE: floor {v['floor']:.3e} > 1%: PERIODIC bar raised from 1% to 2 x floor "
                  f"= {v['bar']:.3e}")
    print("missing/pending (no snapshots):", missing or "none")
    if res:
        OUT.mkdir(parents=True, exist_ok=True)
        (OUT / f"{a.tag}.json").write_text(json.dumps({"runs": res, "missing": missing}))
        make_figure(res, OUT / f"{a.tag}.png")


if __name__ == "__main__":
    main()
