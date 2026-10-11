"""Test I (islands): is the residual tracer variance trapped in a compact, persistent region?

Pre-registration (diary.md, "PRE-REGISTERED test I (islands)", verbatim):

  "5 L8 recordings ... tracer snapshots at 4 frames/period from cycle 80 to 155: rpm 22.5 (slow),
  22.51 (fast; = p8 d1e-2), 22.4 (fast), 30 (base), 29.9 (slow). Run ids isl_l8_rpm<r>_th7.
  Measure: A50(t) = the smallest fraction of the liquid area that holds 50% of the integral of c'^2
  (c' = c - liquid mean), at the snapshots where chi (from the snapshot, solver-consistent estimator)
  is nearest 0.85, 0.90 and 0.95. Also the centroid of that region at the same rocking phase over 5
  consecutive periods around chi = 0.90.
  ISLAND if, at chi = 0.90 and 0.95, A50(slow) <= A50(fast)/2 at 22.5 rpm (slow = 22.5, fast = 22.51
  and 22.4) AND the slow run's A50 region persists: min Jaccard overlap >= 0.5 across the 5 periods at the same
  phase (AMENDED before data; was a centroid shift < 0.1 L). NOT ISLAND if A50(slow)/A50(fast) is in 0.67-1.5 at both chi levels. Else
  inconclusive. 30 vs 29.9 is reported with the same numbers but is not in the verdict.
  Prerequisite: each recording reproduces its original dtmix_0.95 within 2%."

Implementation choices where the text is silent (all also printed / stored in the JSON "choices"):
 1. chi = the native estimator of scripts/chi_scale.py (sigma2_native with the solver's sigma0 and
    f_mean from _solver_reference), identical to chi_series' native chi.
 2. Liquid weight w = cs * f (cell liquid volume; dx^2 is a constant and cancels).  c-bar = sum(w c)/sum(w),
    density = w c'^2, A50 = sum(w of top cells) / sum(w).  The solver's second moment weights interface
    cells by f^2, not f (see chi_scale docstring); w = cs*f is the liquid volume the solver's own
    normaliser V = sum(dv f) uses.  Cells with w = 0 never enter.
 3. A50 uses the smallest set of cells (sorted by density, descending) whose cumulative density first
    reaches 50%; no interpolation inside the last cell.
 4. Persistence: 5 snapshots at times nearest t0 + m*T_per_nd, m = -2..2, t0 = time of the snapshot nearest
    chi = 0.90 (centred on it; shifted inwards if the run ends/starts too early).  T_per_nd from
    postprocess._t_scales.  Centroid = w-weighted mean (x, y) of the top-A50 cells, in L0 units.
    Reported number = MAX PAIRWISE centroid distance over the 5 / L (max distance from the first also stored).
 5. L = 1 in L0 units (domain width; the bag spans the full domain width, measured on a snapshot:
    cs > 0.5 for x/L0 in [0, 0.996]).  geometry.a = 0.25 m = L_bio is the dimensional scale of L0 = 1
    (BioReactor.c comments call it "half-length", but the embedded bag fills the whole L0 = 1 width).
 6. Verdict: "A50(slow) <= A50(fast)/2" must hold against EACH fast run (22.51 and 22.4); NOT ISLAND needs
    the ratio in [0.67, 1.5] against EACH fast run at BOTH chi levels.  The A50 compared is that of the
    snapshot nearest the chi target (not the median of 3).  PENDING = any of the three 22.x runs missing;
    NaN (target never reached) -> INCONCLUSIVE.
 7. 30 vs 29.9: the pre-registration does not say which is "fast"; ratio A50(29.9)/A50(30) is reported.

Usage: uv run python scripts/test_i_islands.py [--runs ID ...]
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import chi_scale as cs_mod                                  # noqa: E402
from scripts.chi_scale import (read_snapshot, locate, sigma2_native,     # noqa: E402
                               _solver_reference, SCRATCH)
from scripts import figstyle as fs                                       # noqa: E402

OUT = ROOT / "experiments" / "multifidelity"
TARGETS = (0.85, 0.90, 0.95)
DEFAULT_RUNS = [f"isl_l8_rpm{r}_th7" for r in ("22.5", "22.51", "22.4", "30", "29.9")]
ORIGINAL = {"isl_l8_rpm22.5_th7": "kmix_l8_rpm22.5", "isl_l8_rpm22.51_th7": "p8_rpm22.5_d1e-2_th7",
            "isl_l8_rpm22.4_th7": "t5_l8_rpm22.4_th7", "isl_l8_rpm30_th7": "kmix_l8_rpm30",
            "isl_l8_rpm29.9_th7": "t5_l8_rpm29.9_th7"}
SLOW, FASTS = "22.5", ("22.51", "22.4")
L_BAG = 1.0            # in L0 units, see choice 5
N_PERIODS = 5
ISLAND_FACTOR, NOT_LO, NOT_HI, JACC_MIN = 0.5, 0.67, 1.5, 0.5
# Persistence (amended 2026-10-11 before any isl data, disclosed in diary): the centroid of a multi-blob A50 set
# is a midpoint that stays put whatever the blobs do (seen on scal_l8 32.5 rpm: two symmetric blobs), so the
# centroid bar had no power. Now: min Jaccard overlap of the A50 masks at the same phase, each period vs the middle
# one, must be >= JACC_MIN. Two unrelated 10%-area sets overlap ~0.05; centroid shift is still reported.


# --------------------------------------------------------------- pure functions
def _density(c, w):
    c = np.asarray(c, float)
    w = np.asarray(w, float)
    cbar = (w * c).sum() / w.sum()
    return w * (c - cbar) ** 2


def a50_mask(c, w):
    """Boolean mask of the top-A50 cell set (same shape as c)."""
    d = _density(c, w)
    w = np.asarray(w, float)
    flat = d.ravel()
    tot = flat.sum()
    mask = np.zeros(flat.shape, bool)
    if not tot > 0:
        return mask.reshape(d.shape)
    order = np.argsort(-flat, kind="stable")
    k = int(np.searchsorted(np.cumsum(flat[order]), 0.5 * tot * (1 - 1e-12))) + 1
    mask[order[:k]] = True
    return mask.reshape(d.shape)


def a50(c, w) -> float:
    """Smallest liquid-area fraction holding 50% of sum w c'^2 (c' about the w-weighted mean).
    NaN if the variance is zero."""
    d = _density(c, w)
    if not d.sum() > 0:
        return math.nan
    m = a50_mask(c, w)
    return float(np.asarray(w, float)[m].sum() / np.asarray(w, float).sum())


def island_rule(slow, fasts, persist_jaccard) -> str:
    """slow = {"0.90": A50, "0.95": A50} or None (run missing); fasts = list of the same (None = missing);
    persist_jaccard = min Jaccard overlap of the slow run's A50 masks across periods (NaN if unknown)."""
    if slow is None or any(f is None for f in fasts):
        return "PENDING"
    levels = ("0.90", "0.95")
    vals = [slow[k] for k in levels] + [f[k] for f in fasts for k in levels]
    if any(v is None or not math.isfinite(v) for v in vals):
        return "INCONCLUSIVE"
    island_a = all(slow[k] <= ISLAND_FACTOR * f[k] for f in fasts for k in levels)
    if island_a and persist_jaccard is not None and math.isfinite(persist_jaccard) \
            and persist_jaccard >= JACC_MIN:
        return "ISLAND"
    if all(NOT_LO <= slow[k] / f[k] <= NOT_HI for f in fasts for k in levels):
        return "NOT ISLAND"
    return "INCONCLUSIVE"


# --------------------------------------------------------------------- data
def rpm_of(run_id):
    return re.search(r"rpm([\d.]+?)_th", run_id).group(1)


def _liquid_w(snap):
    return snap["cs"] * snap["f"]


def analyse_run(run_id, verbose=True):
    loc = locate(run_id)
    if loc["fields_dir"] is None or loc["data_dir"] is None:
        return {"run_id": run_id, "status": "missing"}
    ref = _solver_reference(loc)
    files = sorted(Path(loc["fields_dir"]).glob("snap_*.bin"))
    keep = []                      # (path, t_nd, chi, a50)
    sig0 = ref["sigma0"]
    started = False
    for p in files:
        s = read_snapshot(p, ["f", "cs", "c2"])
        f, cs, c = s["f"], s["cs"], s["c2"]
        if not started:
            if float((cs * f * c).sum()) <= 1e-10 * float((cs * f).sum()):
                continue
            started = True
            if sig0 is None:
                sig0 = sigma2_native(f, cs, c, ref["f_mean"])
        chi = float(np.clip(1.0 - sigma2_native(f, cs, c, ref["f_mean"]) / sig0, 0.0, 1.0))
        keep.append((p, s["t"], chi, a50(c, cs * f)))
        if chi > 0.99:
            break
    if not keep:
        return {"run_id": run_id, "status": "missing", "detail": "no tracer snapshots"}
    t_nd = np.array([k[1] for k in keep])
    chi = np.array([k[2] for k in keep])
    A = np.array([k[3] for k in keep])
    T_bio, T_per = ref["T_bio"], ref["T_per_nd"]
    out = {"run_id": run_id, "status": "ok", "rpm": rpm_of(run_id), "n_snapshots": len(keep),
           "T_per_nd": T_per, "frame_dt_over_T_per": float(np.median(np.diff(t_nd)) / T_per),
           "t_s": ((t_nd - ref["t0"]) * T_bio).tolist(), "chi": chi.tolist(), "a50": A.tolist(),
           "dtmix_0.95_results": (ref["dtmix"] or {}).get("dtmix_0.95"), "targets": {}}
    for tg in TARGETS:
        if chi.max() < tg:
            out["targets"][f"{tg:.2f}"] = {"reached": False, "chi_max": float(chi.max())}
            continue
        k = int(np.argmin(np.abs(chi - tg)))
        near = np.argsort(np.abs(chi - tg))[:3]
        out["targets"][f"{tg:.2f}"] = {
            "reached": True, "index": k, "chi": float(chi[k]), "t_s": out["t_s"][k],
            "t_nd": float(t_nd[k]), "a50": float(A[k]), "file": str(keep[k][0]),
            "a50_3nearest": [float(A[i]) for i in near], "chi_3nearest": [float(chi[i]) for i in near],
            "a50_min": float(A[near].min()), "a50_max": float(A[near].max())}
    out["persistence"] = persistence(keep, t_nd, T_per, out["targets"].get("0.90"))
    out["_keep"] = keep
    return out


def jaccard(a, b) -> float:
    a, b = np.asarray(a, bool), np.asarray(b, bool)
    u = np.logical_or(a, b).sum()
    return float(np.logical_and(a, b).sum() / u) if u else 0.0


def _a50_mask_of(snap_path):
    s = read_snapshot(snap_path, ["f", "cs", "c2"])
    return a50_mask(s["c2"], _liquid_w(s))


def _centroid(snap_path):
    s = read_snapshot(snap_path, ["f", "cs", "c2"])
    w = _liquid_w(s)
    m = a50_mask(s["c2"], w)
    n = s["n"]
    j, i = np.nonzero(m)
    ww = w[j, i]
    x = -0.5 + (i + 0.5) / n
    y = -0.5 + (j + 0.5) / n
    return float((ww * x).sum() / ww.sum()), float((ww * y).sum() / ww.sum())


def persistence(keep, t_nd, T_per, tgt):
    if not tgt or not tgt.get("reached"):
        return {"status": "chi 0.90 not reached"}
    k0 = tgt["index"]
    t0 = t_nd[k0]
    ms = [-2, -1, 0, 1, 2]
    while t0 + ms[0] * T_per < t_nd[0] - 0.5 * T_per and ms[-1] < 8:      # shift inwards
        ms = [m + 1 for m in ms]
    while t0 + ms[-1] * T_per > t_nd[-1] + 0.5 * T_per and ms[0] > -8:
        ms = [m - 1 for m in ms]
    idx = [int(np.argmin(np.abs(t_nd - (t0 + m * T_per)))) for m in ms]
    phase_err = [float((t_nd[i] - (t0 + m * T_per)) / T_per) for i, m in zip(idx, ms)]
    if len(set(idx)) < N_PERIODS or max(abs(e) for e in phase_err) > 0.25:
        return {"status": "not enough snapshots at the same phase", "indices": idx}
    cen = np.array([_centroid(keep[i][0]) for i in idx])
    masks = [_a50_mask_of(keep[i][0]) for i in idx]
    mid = len(masks) // 2
    jac = [jaccard(masks[mid], m) for k, m in enumerate(masks) if k != mid]
    d = np.linalg.norm(cen[:, None] - cen[None], axis=-1)
    return {"status": "ok", "periods_m": ms, "indices": idx, "chi": [float(keep[i][2]) for i in idx],
            "phase_error_over_T": phase_err, "centroids_L0": cen.tolist(),
            "max_pairwise_over_L": float(d.max() / L_BAG),
            "max_from_first_over_L": float(d[0].max() / L_BAG), "L0_units": L_BAG,
            "jaccard_vs_middle": jac, "min_jaccard": float(min(jac))}


def _dtmix_orig(orig):
    for d in (ROOT / "runs" / orig, SCRATCH / orig):
        if (d / "results.json").exists():
            return json.loads((d / "results.json").read_text()).get("dtmix_0.95")
    return None


def prerequisite(res):
    orig = ORIGINAL.get(res["run_id"])
    if orig is None:
        return None
    a, b = res.get("dtmix_0.95_results"), _dtmix_orig(orig)
    if a is None or b is None:
        return {"original": orig, "dtmix_rec": a, "dtmix_orig": b, "within_2pct": None}
    return {"original": orig, "dtmix_rec": a, "dtmix_orig": b, "rel_diff": a / b - 1,
            "within_2pct": bool(abs(a / b - 1) <= 0.02)}


# ------------------------------------------------------------------- figures
def plot_maps(results, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update(fs.rcparams())
    items = []
    for r in results:
        t = r["targets"].get("0.90")
        if t and t.get("reached"):
            s = read_snapshot(t["file"], ["f", "cs", "c2"])
            w = _liquid_w(s)
            cp2 = _density(s["c2"], np.where(w > 0, w, 1e-300)) / np.where(w > 0, w, 1e-300)
            liquid = (s["cs"] > 0.5) & (s["f"] > 0.5)
            items.append((r["rpm"], s, np.where(liquid, cp2, np.nan), a50_mask(s["c2"], w), liquid))
    if not items:
        return False
    groups = {}
    for it in items:
        groups.setdefault(it[0].split(".")[0][:1] + "x", []).append(it)
    n = len(items)
    f, axs = plt.subplots(1, n, figsize=(3.6 * n, 3.2), squeeze=False)
    axs = axs[0]
    k = 0
    for g in groups.values():
        vmax = max(np.nanpercentile(it[2], 99.5) for it in g)
        for rpm, s, field, mask, liquid in g:
            ax = axs[k]
            k += 1
            im = ax.imshow(field, origin="lower", extent=(-0.5, 0.5, -0.5, 0.5), vmin=0, vmax=vmax,
                           cmap="viridis", interpolation="nearest")
            ax.contour(np.linspace(-0.5, 0.5, s["n"]), np.linspace(-0.5, 0.5, s["n"]),
                       mask.astype(float), levels=[0.5], colors="w", linewidths=0.8)
            ax.set_title(f"{rpm} rpm", fontsize=9)
            ax.set_xlabel(r"$x/L_0$")
            ax.set_ylim(-0.2, 0.2)
            ax.set_aspect("equal")
            if k == 1:
                ax.set_ylabel(r"$y/L_0$")
        f.colorbar(im, ax=axs[k - len(g):k], shrink=0.7, pad=0.02, label=r"$c'^2$")
    f.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(f)
    return True


def plot_a50(results, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update(fs.rcparams())
    f, ax = plt.subplots(figsize=(5.2, 3.4))
    markers = ["o", "s", "^", "D", "v", "P"]
    lss = [":", "--", "-.", (0, (3, 1, 1, 1)), (0, (1, 3))]
    for r, mk, ls in zip(results, markers, lss):
        ax.plot(r["chi"], r["a50"], **fs.series_kw(fs.level_colour(8), mk, stat="mean", ours=True,
                                                   ls=ls, ms=3.5, lw=1.0), label=f"{r['rpm']} rpm")
    ax.set_xlabel(r"$\chi$ (-)")
    ax.set_ylabel(r"$A_{50}$ (-)")
    ax.set_yscale("log")
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=8,
              title="L8; line style/marker = rpm", title_fontsize=8)
    f.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(f)


# ---------------------------------------------------------------------- main
def _a50_at(res, tg):
    if res is None or res["status"] != "ok":
        return None
    t = res["targets"].get(tg)
    return t["a50"] if t and t.get("reached") else math.nan


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="*", default=DEFAULT_RUNS)
    args = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    results, missing = [], []
    for rid in args.runs:
        r = analyse_run(rid)
        if r["status"] != "ok":
            missing.append(rid)
            print(f"[{rid}] MISSING (no tracer snapshots)")
            continue
        results.append(r)
        print(f"\n[{rid}] {r['n_snapshots']} snapshots, frame dt = {r['frame_dt_over_T_per']:.4f} T_per "
              f"(T_per_nd = {r['T_per_nd']:.4f})")
        for tg, t in r["targets"].items():
            if not t["reached"]:
                print(f"  chi target {tg}: NOT REACHED (max chi {t['chi_max']:.4f})")
                continue
            print(f"  chi target {tg}: hit chi = {t['chi']:.4f} at t = {t['t_s']:.2f} s "
                  f"(snapshot {t['index']}): A50 = {t['a50']:.4f}; 3 nearest A50 min/max = "
                  f"{t['a50_min']:.4f} / {t['a50_max']:.4f}")
        p = r["persistence"]
        if p.get("status") == "ok":
            print(f"  persistence (periods {p['periods_m']}, chi {np.round(p['chi'], 3).tolist()}): "
                  f"max pairwise centroid shift / L = {p['max_pairwise_over_L']:.4f} "
                  f"(from first: {p['max_from_first_over_L']:.4f}); Jaccard vs middle = "
                  f"{np.round(p['jaccard_vs_middle'], 3).tolist()} (min {p['min_jaccard']:.3f}); max |phase err| = "
                  f"{max(abs(e) for e in p['phase_error_over_T']):.3f} T")
        else:
            print(f"  persistence: {p['status']}")
        pr = prerequisite(r)
        r["prerequisite"] = pr
        if pr:
            if pr["within_2pct"] is None:
                print(f"  prerequisite vs {pr['original']}: original/recording dtmix_0.95 unavailable")
            else:
                print(f"  prerequisite vs {pr['original']}: dtmix_0.95 {pr['dtmix_rec']:.3f} vs "
                      f"{pr['dtmix_orig']:.3f} ({100 * pr['rel_diff']:+.2f}%) -> "
                      f"{'OK' if pr['within_2pct'] else 'FAIL (reassign by own chi(t))'}")
    by = {r["rpm"]: r for r in results}
    miss_needed = [f"isl_l8_rpm{r}_th7" for r in (SLOW,) + FASTS if r not in by]
    levels = ("0.90", "0.95")

    def lv(rpm):
        r = by.get(rpm)
        return None if r is None else {k: _a50_at(r, k) for k in levels}
    slow = lv(SLOW)
    persist = by[SLOW]["persistence"].get("min_jaccard", math.nan) if SLOW in by else math.nan
    verdict = island_rule(slow, [lv(f) for f in FASTS], persist)
    print("\n=== VERDICT (22.5 slow vs 22.51, 22.4 fast; chi 0.90 and 0.95) ===")
    print(f"{verdict}" + (f"   missing runs: {miss_needed}" if miss_needed else ""))
    ratios = {}
    for f in FASTS:
        if slow and lv(f):
            ratios[f] = {k: slow[k] / lv(f)[k] for k in levels}
            print(f"  A50(22.5)/A50({f}): " + ", ".join(f"chi {k}: {v:.3f}" for k, v in ratios[f].items()))
    pair30 = None
    if "30" in by and "29.9" in by:
        pair30 = {k: _a50_at(by["29.9"], k) / _a50_at(by["30"], k) for k in levels}
        print("  (report only) A50(29.9)/A50(30): " + ", ".join(f"chi {k}: {v:.3f}" for k, v in pair30.items()))
    if missing:
        print(f"\nPENDING / missing runs: {missing}")
    figs = {}
    if results:
        figs["maps"] = plot_maps(results, OUT / "test_i_islands.png")
        plot_a50(results, OUT / "test_i_a50.png")
        print(f"\nfigures: {OUT / 'test_i_islands.png'} ({'written' if figs['maps'] else 'no chi 0.90 snapshot'}), "
              f"{OUT / 'test_i_a50.png'}")
    summ = {"verdict": verdict, "missing_runs": missing, "missing_needed_for_verdict": miss_needed,
            "ratios_slow_over_fast": ratios, "a50_29.9_over_30": pair30,
            "L_bag_L0_units": L_BAG, "rule": {"island_factor": ISLAND_FACTOR,
                                               "not_island_ratio": [NOT_LO, NOT_HI], "min_jaccard": JACC_MIN},
            "choices": [l.strip() for l in __doc__.split("Implementation choices")[1].split("Usage")[0].splitlines()
                        if l.strip()],
            "runs": [{k: v for k, v in r.items() if k != "_keep"} for r in results]}
    (OUT / "test_i_islands.json").write_text(json.dumps(cs_mod._json_safe(summ), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
