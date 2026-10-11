"""Degree of mixing chi(t) at the solver's grid scale and at fixed physical scales.

Purpose (spec, diary.md, quoted verbatim)
-----------------------------------------
PRE-REGISTERED test L:

  "- PRE-REGISTERED test L (scale artefact vs not), compute ~1,000 core-h: L8
  and L9, theta 7, 32.5 and 37.5 rpm, lean binary f1c11e0, c snapshots
  13/period over cycles 80-180 (scal_l<L>_rpm<r>_th7). chi_l(t) from c2 by
  Kim's PROJECT coarse-graining onto fixed boxes l = L0/16, L0/32, L0/64
  (liquid cells, volume-fraction weighted), and at the native grid scale
  (check: native must reproduce results.json dtmix within 2%).
  Reading, per rpm and chi in {0.50, 0.75, 0.95}, d_grid = |log dtmix_L9 - log
  dtmix_L8| at native scale, d_l the same at scale l:
  SCALE ARTEFACT if d_l <= d_grid / 3 for l = L0/32 at >= 5 of the 6 (rpm,
  chi) cases; NOT an artefact if d_l >= 2 d_grid / 3 at >= 5 of 6; otherwise
  mixed (reported as such)."

Test R threshold REVISED:

  "- Test R threshold REVISED before any rec_l8 result exists (disclosed):
  H_rep holds if the median |log rerun - log original| of dtmix_0.95 over the
  8 reruns >= 0.10 (1/2 of the median +-0.1-rpm jump, ~0.20); reruns are
  reproducible if it is <= 0.05. Old threshold (0.02 at >= 2 of 8) dropped: it
  passes trivially for small effects."

What the solver's own chi is (read from src/BioReactor.c normcal and
scripts/postprocess.py::_compute_mixing_metrics)
------------------------------------------------
* Tracer: c2 (top half), plane "c2" of the snapshot; tr_oxy.dat cols 8/9.
* Per cell: a = c2 * f   (f = liquid volume fraction), weight dv = cs * dx^2
  (cs = embedded fluid fraction; cs = 1 in fluid cells, 0 in solid).
  sum = Sum dv a, sum2 = Sum dv a^2 (statsf2).  NOTE sum2 is Sum dv (c f)^2,
  i.e. f enters SQUARED in the second moment, not Sum dv f c^2.
* Normalisation: c_mean = sum / f_mean, c_mean2 = sum2 / f_mean with f_mean the
  TIME-MEAN over the whole run of f_liq_sum (vol_frac_interf.dat col 2),
  = Sum dv f.  sigma2 = c_mean2 - c_mean^2.
* sigma0 = sigma2 at the first tr_oxy.dat row with non-zero tracer (the
  injection, logged at the t_out cadence); chi = clip(1 - sigma2/sigma0, 0, 1);
  dtmix = (t of FIRST logged row with chi >= thr  -  t of that first row) * T_bio.
  No interpolation in time.

Kim's postprocessing (dev/postprocessing/bio_mix_elvis.m) differs
------------------------------------------------------------------
* Cells/boxes counted as liquid when the (projected) alpha > 0, then every
  such cell/box is weighted EQUALLY (dA uniform); no f weighting, no solid
  fraction weighting; tracer is the raw (block-mean) value, not c*f.
* sigma0 = first sample at/after t_mix (in1 == 1).
* With PROJECT = 1 the raw tracer is block-averaged over all cells of a box
  (including gas and solid cells), then the above.
So the solver's chi and Kim's chi are not the same estimator even at the
native scale; both are implemented here ("solver"/native, "vol", "kim").

Estimators implemented
----------------------
native   solver's definition above (sigma2_native, uses c*f, f squared).
vol      PRIMARY coarse-graining: conservative liquid-volume-weighted box mean
         C_b = Sum cs f c / Sum cs f over the box, box weight V_b = Sum cs f;
         sigma2 = Sum V_b C_b^2 / V - (Sum V_b C_b / V)^2, V = Sum V_b.  Boxes
         with V_b = 0 contribute nothing (excluded).  At block = 1 this equals
         the solver's sigma2 exactly where f in {0,1} and differs only in the
         interface cells (f in (0,1)): solver has f^2 c^2, this has f c^2.
kim      SECONDARY: Kim's raw block mean of c over all cells of the box;
         boxes kept if their mean f > 0; unweighted; sigma2 = <C^2> - <C>^2.

Choices where the spec was silent
---------------------------------
1. Coarse-graining of the solver's quantity: the spec says "volume-weighted,
   follow the solver".  Weighting the solver's a = c f by box VOLUME (instead of
   by liquid volume) would put a spurious variance at the free surface that
   grows with box size (a half-gas box has a_b = c/2).  So "vol" averages c
   (not c f) with weight cs f.  Consequence: block = 1 matches the solver only
   up to the interface-cell treatment (see above); the gap is measured on real
   data by compare_interface_treatment() and reported by the CLI.
2. Native chi from snapshots uses the SOLVER's sigma0 (results/tr_oxy sigma2_max)
   and the solver's f_mean (time-mean of vol_frac_interf.dat) when those files
   exist; else sigma0 of the first snapshot with tracer and the snapshot's own
   Sum dv f.  Coarse sigma0_l always comes from the first snapshot that has
   tracer, because the coarse-grained variance cannot be had from the .dat
   files.  If that snapshot is later than the release row, its native chi is
   recorded as "release_offset_chi" so the bias is visible.
3. Time origin for dtmix at every scale = t of the first non-zero tracer row of
   tr_oxy.dat (the solver's own origin); fallback t_checkpoint + n_mix_cycles *
   T_per_nd.  A point (t0, chi = 0) is prepended to the series when the first
   snapshot is after the origin.
4. Crossing time: linear interpolation of chi in time between the two snapshots
   bracketing the threshold.  Also stored: first-sample crossing
   (dtmix_native_first) to separate cadence effects from definition effects.
5. Reading stops once every chi series (native + all scales, both estimators)
   has exceeded 0.99 (nothing later can move a dtmix <= 0.95).
6. Test L uses the SNAPSHOT-derived native dtmix for d_grid (same estimator and
   cadence as d_l); the results.json-based d_grid is stored alongside as a
   sensitivity. The primary l estimator for the verdict is "vol" at L0/32; the
   verdict with Kim's raw mean is stored as a secondary.
7. Verdict with a missing run: "undetermined (runs missing)".  A case whose
   threshold is never reached by one of the runs is neither counted nor
   "missing".  The amended rule (below) is implemented as: counted = cases with
   d_grid >= 0.20 (snapshot-derived native); ARTEFACT if counted - n_artefact
   <= 1; NOT if counted - n_not <= 1; else mixed; < 2 counted -> inconclusive.
   (If both "all but one" conditions held, ARTEFACT is reported first; they can
   only both hold with d_l ~ 0 and d_l >= 2/3 d_grid together, i.e. d_grid ~ 0,
   excluded by the 0.20 floor.)
8. Test R (amended): original = runs/kmix_l8_rpm<r> if present (grid rpm), else
   runs/t5_l8_rpm<r>_th7; per-run |dlog| for the three dtmix, exactly-zero flag.
9. Test P: the offset-0.1 point uses the +0.1 run (t5 at base + 0.1) so that
   its sign matches the p8 offsets; the -0.1 run is stored, not used in the
   reading.  "fall by >= 10x per 2 decades" is implemented as
   med(0.1)/med(1e-3) >= 10 and med(1e-2)/med(1e-6) >= 100 (zero denominator
   satisfies it).  Medians are over the two bases and require both.  Exact
   zeros are drawn at 1e-8 on the log axis.  Base runs: kmix_l8_rpm22.5 and
   kmix_l8_rpm30.
10. Figure (test L): ours dotted (figstyle); native = thin line, l = L0/32 (vol) = thick
   line; L8/L9 colours from figstyle; only native and L0/32 are drawn (other
   scales are in the JSON); legend outside axes at right; one panel per rpm.
   Test P figure: the two bases get two non-level colours (blue/vermillion
   Okabe-Ito) because figstyle reserves level hues.

Spec amendment (coordinator, 2026-10-10, diary.md "POWER of the pending tests"),
quoted:

  "Test R (exact reruns) ... R is kept only as a determinism check (does the T5
  jump include run-to-run nondeterminism?). It is no longer the test of H_rep."
  "NEW test P (replaces R as the H_rep test), PRE-REGISTERED: perturbation-size
  scan at L8, theta 7, base 22.5 and 30 rpm (T5 +-0.1 jumps 0.20-0.35 in log),
  offsets +1e-2, +1e-3, +1e-6 rpm (6 runs, scripts/submit_l8_perturb.py; with
  the T5 +-0.1 runs that gives 4 decades).
  Model: |Delta log dtmix_0.95| vs |Delta rpm|. Smooth-but-steep structure
  predicts slope ~1 in log-log ... Sensitive dependence (effectively noise)
  predicts no decay ...
  Reading: NOISE if the median |Delta| over the two bases at 1e-6 rpm >= 0.10;
  SMOOTH if at 1e-3 rpm it is <= 0.02 AND the jumps fall by >= 10x per 2
  decades; otherwise mixed."
  "Test L power, amended: ... New rule: only (rpm, chi) cases with d_grid >=
  0.20 in the NEW scal runs count (expected 3-4 cases: 32.5 rpm chi
  0.50/0.75/0.95, 37.5 rpm chi 0.95). SCALE ARTEFACT if d_l <= d_grid/3 at l =
  L0/32 in all but at most one counted case; NOT if d_l >= 2 d_grid/3 in all
  but at most one; else mixed. Fewer than 2 counted cases -> test L
  inconclusive."
  Test G (PRE-REGISTERED, diary.md line "PRE-REGISTERED test G"), verbatim:
  "- PRE-REGISTERED test G (scalar gradient scale, zero extra compute; data from
  scal_l8/scal_l9 and rec_l8 snapshots): replaces my hand-waved "sheets of pure
  dye" picture with a measurement. Scalar microscale lambda_c(t) = sqrt(<c'^2> /
  <|grad c|^2>) over liquid cells (central differences on the uniform snapshot
  grid), after release.
  GRID-CONTROLLED if, during the decay (chi 0.5 -> 0.95), lambda_c/dx is O(1-3)
  at both L8 and L9 AND lambda_c(L8)/lambda_c(L9) is within 1.5-2.5 (i.e.
  lambda_c scales with dx). NOT grid-controlled if lambda_c is the same physical
  length at both levels (ratio 0.8-1.25) and >= 5 dx.
  Power: the two outcomes differ by a factor 2 in the L8/L9 ratio, measured from
  thousands of cells per snapshot; snapshot-to-snapshot scatter is reported and
  must be below the 0.3 margin.
  Also reported: the fraction of liquid cells with 0.1 < c < 0.9 vs time (direct
  evidence for or against sharp striations)."
  Test G choices: liquid mask f >= 0.999 and cs >= 0.999, gradient needs the cell
  and its 4 neighbours liquid and not on the domain edge; c' about the mean of
  those cells; "O(1-3)" = 1 <= lambda_c/dx <= 3; window = snapshots whose native
  (solver-definition) chi is in [0.5, 0.95]; window medians; dx[m] = 0.25/n;
  overall verdict = per-rpm verdicts if they agree, else "neither/mixed across rpm";
  scatter = std/mean of lambda_c in the window, combined in quadrature.

  SUPERSEDED original test L rule (kept for the record): "SCALE ARTEFACT if d_l
  <= d_grid / 3 for l = L0/32 at >= 5 of the 6 (rpm, chi) cases; NOT an artefact
  if d_l >= 2 d_grid / 3 at >= 5 of 6; otherwise mixed".  SUPERSEDED original
  test R text (H_rep threshold): "H_rep holds if the median |log rerun - log
  original| of dtmix_0.95 over the 8 reruns >= 0.10 ...; reruns are reproducible
  if it is <= 0.05."

  Test R as implemented reports |log rerun - log original| per run for
  dtmix_0.50/0.75/0.95 and whether each is exactly 0; the H_rep verdict (and the
  "Test R threshold REVISED" text quoted at the top) is superseded.  P also
  reports dtmix_0.50/0.75 and kLa_1T_10/25/50 (not criteria).

Usage
-----
  uv run python scripts/chi_scale.py L
  uv run python scripts/chi_scale.py R
  uv run python scripts/chi_scale.py P
  uv run python scripts/chi_scale.py G
  uv run python scripts/chi_scale.py quant
  uv run python scripts/chi_scale.py check RUN_ID        # native vs results.json
"""
from __future__ import annotations

import argparse
import json
import math
import struct
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

PLANE_NAMES = ["ux", "uy", "omega", "f", "cs", "c", "c1", "c2", "c3", "oxy"]
SCALE_DIVISORS = (16, 32, 64)           # l = L0 / divisor
THRESHOLDS = (0.50, 0.75, 0.95)
SCRATCH = Path("/oscar/scratch/eaguerov/mpi_runs")
OUT_DIR = ROOT / "experiments" / "multifidelity"
_HEADER = struct.Struct("<iid")


# ----------------------------------------------------------------- reader
def read_snapshot(path, planes=None) -> dict:
    """Read one solver snapshot (BioReactor.c write_snapshot).

    Layout: int n, int np, double t_nd, then np planes of n*n float32, each
    row-major with j (y, row 0 = bottom) slow and i (x) fast.  `planes` is an
    optional list of plane names; only those are read (memory-mapped).
    Returns {"n", "np", "t", <name>: float64 array (n, n) ...}.
    """
    path = Path(path)
    with open(path, "rb") as fh:
        n, npl, t = _HEADER.unpack(fh.read(_HEADER.size))
    names = PLANE_NAMES[:npl]
    out = {"n": n, "np": npl, "t": t}
    want = names if planes is None else list(planes)
    for name in want:
        k = names.index(name)
        arr = np.memmap(path, dtype="<f4", mode="r",
                        offset=_HEADER.size + k * n * n * 4, shape=(n, n))
        out[name] = np.array(arr, dtype=np.float64)
    return out


def iter_snapshots(fields_dir, planes=None):
    """Yield snapshots of a fields/ directory in time order (file index order,
    which is time order; verified by sorting on the header time)."""
    files = sorted(Path(fields_dir).glob("snap_*.bin"))
    for p in files:
        yield read_snapshot(p, planes)


# ------------------------------------------------------------- estimators
def sigma2_native(f, cs, c, f_mean=None) -> float:
    """Solver's variance: (Sum dv (c f)^2)/V - (Sum dv c f / V)^2, dv = cs dx^2.
    V = f_mean (time-mean of f_liq_sum) if given, else Sum dv f."""
    n = f.shape[0]
    dv = cs / (n * n)
    a = c * f
    s1 = float((dv * a).sum())
    s2 = float((dv * a * a).sum())
    v = float((dv * f).sum()) if f_mean is None else float(f_mean)
    return s2 / v - (s1 / v) ** 2


def coarse_grain(f, cs, c, block):
    """Liquid-volume-weighted box means.  Returns (V_b, C_b), each shape
    (n/block, n/block): V_b = Sum cs f over the box, C_b = Sum cs f c / V_b
    (0 where V_b = 0).  Conserves Sum V_b C_b = Sum cs f c."""
    n = f.shape[0]
    if n % block:
        raise ValueError(f"block {block} does not divide n = {n}")
    m = n // block
    w = (cs * f).reshape(m, block, m, block)
    Vb = w.sum(axis=(1, 3))
    num = (w * c.reshape(m, block, m, block)).sum(axis=(1, 3))
    Cb = np.divide(num, Vb, out=np.zeros_like(num), where=Vb > 0)
    return Vb, Cb


def sigma2_coarse(f, cs, c, block, mode="vol") -> float:
    """Variance of the tracer at box scale `block` cells.  mode "vol" =
    primary (liquid-volume weighted); "kim" = raw block mean, unweighted,
    boxes with mean f > 0 (bio_mix_elvis.m PROJECT)."""
    if mode == "vol":
        Vb, Cb = coarse_grain(f, cs, c, block)
        v = Vb.sum()
        m1 = (Vb * Cb).sum() / v
        m2 = (Vb * Cb * Cb).sum() / v
        return float(m2 - m1 * m1)
    if mode == "kim":
        n = f.shape[0]
        m = n // block
        Cb = c.reshape(m, block, m, block).mean(axis=(1, 3))
        Fb = f.reshape(m, block, m, block).mean(axis=(1, 3))
        sel = Fb > 0
        x = Cb[sel]
        return float((x * x).mean() - x.mean() ** 2)
    raise ValueError(mode)


def crossing_time(t, chi, thr, first_sample=False) -> float:
    """Time at which chi first reaches thr; linear interpolation between the
    bracketing samples (or the first sample >= thr if first_sample).  NaN if
    never reached."""
    t = np.asarray(t, float)
    chi = np.asarray(chi, float)
    hit = np.nonzero(chi >= thr)[0]
    if hit.size == 0:
        return math.nan
    k = int(hit[0])
    if first_sample or k == 0:
        return float(t[k])
    t0, t1, c0, c1 = t[k - 1], t[k], chi[k - 1], chi[k]
    return float(t0 + (thr - c0) * (t1 - t0) / (c1 - c0))


# -------------------------------------------------------------- run access
def locate(run_id: str) -> dict:
    """Where a run's files are: fields_dir (None if no snapshots anywhere),
    data_dir (holds tr_oxy.dat/params.json), results.json path."""
    cands = [ROOT / "runs" / run_id, SCRATCH / run_id]
    fields = None
    for d in cands:
        fd = d / "fields"
        if fd.is_dir() and any(fd.glob("snap_*.bin")):
            fields = fd
            break
    data = None
    for d in ([fields.parent] if fields else []) + cands:
        if (d / "params.json").exists() and (d / "tr_oxy.dat").exists():
            data = d
            break
    if data is None:
        for d in cands:
            if (d / "params.json").exists():
                data = d
                break
    res = ROOT / "runs" / run_id / "results.json"
    return {"fields_dir": fields, "data_dir": data,
            "results": res if res.exists() else None}


def _solver_reference(loc: dict) -> dict:
    """Everything the solver's own pipeline says about the run (via
    scripts/postprocess.py, no re-implementation): params, T_bio, sigma0,
    f_mean, release time, and dtmix from results.json (or recomputed)."""
    from scripts import postprocess as pp
    d = loc["data_dir"]
    params = json.loads((d / "params.json").read_text())
    T_bio, T_per_nd = pp._t_scales(params)
    ref = {"params": params, "T_bio": T_bio, "T_per_nd": T_per_nd,
           "sigma0": None, "f_mean": None, "t0": None, "dtmix": {}}
    if (d / "tr_oxy.dat").exists() and (d / "vol_frac_interf.dat").exists():
        try:
            t, c_sum, _ = pp._joined_cols(d, "tr_oxy.dat",
                                          [pp._COL_C2_LIQ_SUM, pp._COL_C2_LIQ_SUM2])
            _, f_liq = pp._joined_cols(d, "vol_frac_interf.dat", [pp._COL_F_LIQ])
            ref["f_mean"] = float(f_liq.mean())
            nz = np.where(c_sum > 1e-10 * ref["f_mean"])[0]
            if nz.size:
                ref["t0"] = float(t[nz[0]])
            m = pp._compute_mixing_metrics(d, params)
            ref["sigma0"] = m.get("sigma2_max")
            ref["dtmix_recomputed"] = {k: v for k, v in m.items() if k.startswith("dtmix")}
        except (FileNotFoundError, ValueError):
            pass
    if loc["results"]:
        r = json.loads(loc["results"].read_text())
        ref["dtmix"] = {k: r.get(k) for k in ("dtmix_0.50", "dtmix_0.75", "dtmix_0.95")}
        if ref["sigma0"] is None:
            ref["sigma0"] = r.get("sigma2_max")
    elif "dtmix_recomputed" in ref:
        ref["dtmix"] = {k: v for k, v in ref["dtmix_recomputed"].items()}
    if ref["t0"] is None:
        ref["t0"] = float(params.get("t_checkpoint") or 0.0) + \
            T_per_nd * float(params.get("n_mix_cycles", 80))
        ref["t0_fallback"] = True
    return ref


# ------------------------------------------------------------ chi of a run
def chi_series(run_id: str, divisors=SCALE_DIVISORS, thresholds=THRESHOLDS,
               verbose=False) -> dict:
    """chi(t) at the native scale and at l = L0/d for each d, from snapshots.

    Returns a dict (status "ok") with t_s (seconds after release), chi
    {"native": arr, "vol_<d>": arr, "kim_<d>": arr}, dtmix {key: {thr: s}},
    solver reference values and diagnostics; or {"status": "missing", ...}.
    """
    loc = locate(run_id)
    if loc["fields_dir"] is None or loc["data_dir"] is None:
        return {"run_id": run_id, "status": "missing",
                "detail": "no fields/snap_*.bin in runs/ or mpi_runs/"}
    ref = _solver_reference(loc)
    ts, s2 = [], {}
    sig0_l, first_snap = {}, None
    for snap in iter_snapshots(loc["fields_dir"], ["f", "cs", "c2"]):
        f, cs, c = snap["f"], snap["cs"], snap["c2"]
        n = snap["n"]
        if first_snap is None:
            if float((cs * f * c).sum()) <= 1e-10 * float((cs * f).sum()):
                continue                    # before injection
            first_snap = snap["t"]
        ts.append(snap["t"])
        row = {"native": sigma2_native(f, cs, c, ref["f_mean"])}
        for d in divisors:
            for mode in ("vol", "kim"):
                row[f"{mode}_{d}"] = sigma2_coarse(f, cs, c, n // d, mode)
        for k, v in row.items():
            s2.setdefault(k, []).append(v)
            sig0_l.setdefault(k, v)         # first snapshot with tracer
        if all(sig0_l[k] > 0 and (1 - s2[k][-1] / sig0_l[k]) > 0.99 for k in s2):
            break
    if not ts:
        return {"run_id": run_id, "status": "missing", "detail": "snapshots without tracer"}
    ts = np.array(ts)
    s2 = {k: np.array(v) for k, v in s2.items()}
    sigma0 = {k: (ref["sigma0"] if k == "native" and ref["sigma0"] else sig0_l[k])
              for k in s2}
    chi = {k: np.clip(1.0 - s2[k] / sigma0[k], 0.0, 1.0) for k in s2}
    t_rel = (ts - ref["t0"]) * ref["T_bio"]
    if t_rel[0] > 0:                         # anchor chi = 0 at release
        t_anchor = np.concatenate([[0.0], t_rel])
        chi_a = {k: np.concatenate([[0.0], v]) for k, v in chi.items()}
    else:
        t_anchor, chi_a = t_rel, chi
    dtmix = {k: {f"{thr:.2f}": crossing_time(t_anchor, chi_a[k], thr) for thr in thresholds}
             for k in chi_a}
    dtmix_first = {f"{thr:.2f}": crossing_time(t_anchor, chi_a["native"], thr, True)
                   for thr in thresholds}
    return {"run_id": run_id, "status": "ok", "n_snapshots": int(ts.size),
            "t_s": t_anchor, "chi": chi_a, "dtmix": dtmix,
            "dtmix_native_first": dtmix_first,
            "solver": {"dtmix": ref["dtmix"], "sigma0": ref["sigma0"],
                       "f_mean": ref["f_mean"], "t0_nd": ref["t0"],
                       "t0_fallback": bool(ref.get("t0_fallback"))},
            "sigma0_used": sigma0,
            "release_offset_s": float((first_snap - ref["t0"]) * ref["T_bio"]),
            "release_offset_chi": float(chi["native"][0]),
            "fields_dir": str(loc["fields_dir"]),
            "t_end_s": float(t_anchor[-1])}


def compare_to_results(series: dict) -> dict:
    """Native dtmix from snapshots vs the solver's (results.json) per threshold."""
    out = {}
    for thr in ("0.50", "0.75", "0.95"):
        ref = series["solver"]["dtmix"].get(f"dtmix_{thr}")
        mine = series["dtmix"]["native"][thr]
        first = series["dtmix_native_first"][thr]
        if ref is None or not np.isfinite(ref) or not np.isfinite(mine):
            out[thr] = {"solver": ref, "snapshots": None if not np.isfinite(mine) else mine,
                        "rel_diff": None, "within_2pct": None}
            continue
        rel = mine / ref - 1.0
        out[thr] = {"solver": float(ref), "snapshots_interp": mine,
                    "snapshots_first_sample": first,
                    "rel_diff": rel, "rel_diff_first_sample": first / ref - 1.0,
                    "within_2pct": bool(abs(rel) <= 0.02)}
    return out


def compare_interface_treatment(run_id: str, n_snap: int = 5) -> dict:
    """On real snapshots: solver sigma2 (c f squared) vs "vol" estimator at
    block 1, relative gap -- the interface-cell difference of choice 1."""
    loc = locate(run_id)
    if loc["fields_dir"] is None:
        return {"status": "missing"}
    files = sorted(loc["fields_dir"].glob("snap_*.bin"))
    pick = files[:: max(1, len(files) // n_snap)][:n_snap]
    gaps = []
    for p in pick:
        s = read_snapshot(p, ["f", "cs", "c2"])
        a = sigma2_native(s["f"], s["cs"], s["c2"])
        b = sigma2_coarse(s["f"], s["cs"], s["c2"], 1, "vol")
        if a > 0:
            gaps.append(b / a - 1.0)
    return {"status": "ok", "rel_gap_vol_vs_solver_block1": gaps}


# --------------------------------------------------------------- test L
L_RUNS = {(8, 32.5): "scal_l8_rpm32.5_th7", (9, 32.5): "scal_l9_rpm32.5_th7",
          (8, 37.5): "scal_l8_rpm37.5_th7", (9, 37.5): "scal_l9_rpm37.5_th7"}


def _dlog(a, b):
    if a is None or b is None or not (np.isfinite(a) and np.isfinite(b)) or a <= 0 or b <= 0:
        return None
    return abs(math.log(b) - math.log(a))


D_GRID_MIN = 0.20      # amended rule: only cases with d_grid >= 0.20 are counted


def test_l_verdict(series: dict, scale_key: str = "vol_32") -> dict:
    """series[(level, rpm)] -> chi_series dict. Applies the pre-registered rule."""
    cases = []
    for rpm in (32.5, 37.5):
        s8, s9 = series[(8, rpm)], series[(9, rpm)]
        for thr in ("0.50", "0.75", "0.95"):
            case = {"rpm": rpm, "chi": thr, "d_grid": None, "d_l": None, "status": "missing"}
            if s8["status"] == "ok" and s9["status"] == "ok":
                case["d_grid"] = _dlog(s8["dtmix"]["native"][thr], s9["dtmix"]["native"][thr])
                case["d_l"] = _dlog(s8["dtmix"][scale_key][thr], s9["dtmix"][scale_key][thr])
                r8 = s8["solver"]["dtmix"].get(f"dtmix_{thr}")
                r9 = s9["solver"]["dtmix"].get(f"dtmix_{thr}")
                case["d_grid_results_json"] = _dlog(r8, r9)
                case["d_l_over_d_grid"] = (case["d_l"] / case["d_grid"]
                                           if case["d_grid"] and case["d_l"] is not None else None)
                if case["d_grid"] is not None and case["d_l"] is not None:
                    case["status"] = "ok"
                    case["artefact"] = bool(case["d_l"] <= case["d_grid"] / 3)
                    case["not_artefact"] = bool(case["d_l"] >= 2 * case["d_grid"] / 3)
                else:
                    case["status"] = "unreached"      # listed, not counted
                    un = []
                    for lvl, sr in ((8, s8), (9, s9)):
                        for k in ("native", scale_key):
                            v = sr["dtmix"][k][thr]
                            if v is None or not math.isfinite(v):
                                arr = sr.get("chi", {}).get(k)
                                un.append({"level": lvl, "estimator": k,
                                           "max_chi_reached": None if arr is None else float(np.max(arr))})
                    case["unreached"] = un
            cases.append(case)
    missing = [c for c in cases if c["status"] == "missing"]
    ok = [c for c in cases if c["status"] == "ok"]
    for c in ok:
        c["counted"] = bool(c["d_grid"] >= D_GRID_MIN)
    counted = [c for c in ok if c["counted"]]
    n_art = sum(c["artefact"] for c in counted)
    n_not = sum(c["not_artefact"] for c in counted)
    if missing:
        verdict = "undetermined (runs missing)"
    elif len(counted) < 2:
        verdict = "inconclusive (fewer than 2 counted cases)"
    elif len(counted) - n_art <= 1:
        verdict = "SCALE ARTEFACT"
    elif len(counted) - n_not <= 1:
        verdict = "NOT an artefact"
    else:
        verdict = "mixed"
    return {"scale": scale_key, "d_grid_min": D_GRID_MIN, "n_counted": len(counted),
            "n_artefact": n_art, "n_not_artefact": n_not,
            "n_missing": len(missing), "verdict": verdict, "cases": cases}


def _json_safe(o):
    if isinstance(o, dict):
        return {str(k): _json_safe(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_json_safe(v) for v in o]
    if isinstance(o, np.ndarray):
        return [None if not np.isfinite(x) else float(x) for x in o]
    if isinstance(o, (np.floating, float)):
        return None if not math.isfinite(o) else float(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def plot_test_l(series: dict, path: Path) -> bool:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scripts import figstyle as fs
    plt.rcParams.update(fs.rcparams())
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
    drew = False
    for ax, rpm in zip(axes, (32.5, 37.5)):
        for lvl in (8, 9):
            s = series[(lvl, rpm)]
            if s["status"] != "ok":
                continue
            col = fs.level_colour(lvl)
            ax.plot(s["t_s"], s["chi"]["native"], color=col, ls=":", lw=1.0,
                    label=f"L{lvl}, grid scale")
            ax.plot(s["t_s"], s["chi"]["vol_32"], color=col, ls=":", lw=2.6,
                    label=f"L{lvl}, $l=L_0/32$")
            drew = True
        ax.set_title(f"{rpm:g} rpm, 7$^\\circ$")
        ax.set_xlabel("time after release (s)")
        ax.grid(**fs.GRID_KW)
    axes[0].set_ylabel("degree of mixing $\\chi$ (-)")
    if drew:
        h, l = axes[1].get_legend_handles_labels()
        axes[1].legend(h, l, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    if drew:
        fig.tight_layout()
        fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return drew


def cmd_L(args) -> int:
    series = {key: chi_series(rid) for key, rid in L_RUNS.items()}
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    res = {"runs": {}, "verdict_primary": None, "verdict_kim_secondary": None}
    for (lvl, rpm), s in series.items():
        entry = {k: v for k, v in s.items() if k not in ("t_s", "chi")}
        if s["status"] == "ok":
            entry["native_vs_results_json"] = compare_to_results(s)
            entry["interface_treatment"] = compare_interface_treatment(s["run_id"])
        res["runs"][s["run_id"]] = entry
    res["verdict_primary"] = test_l_verdict(series, "vol_32")
    res["verdict_kim_secondary"] = test_l_verdict(series, "kim_32")
    res["other_scales"] = {k: test_l_verdict(series, k)["verdict"]
                           for k in ("vol_16", "vol_64", "kim_16", "kim_64")}
    (OUT_DIR / "test_l_scale.json").write_text(json.dumps(_json_safe(res), indent=1))
    plot_test_l(series, OUT_DIR / "test_l_scale.png")
    for key, s in series.items():
        print(key, s["status"], s.get("detail", ""))
    print("verdict (vol, L0/32):", res["verdict_primary"]["verdict"])
    return 0


# --------------------------------------------------------------- test R
R_RPMS = ["22.4", "22.5", "22.6", "25", "29.9", "30", "30.1", "32.5"]
DT_KEYS = ("dtmix_0.50", "dtmix_0.75", "dtmix_0.95")
KLA_KEYS = ("kLa_1T_10", "kLa_1T_25", "kLa_1T_50")


def _res(run_id: str):
    p = ROOT / "runs" / run_id / "results.json"
    return json.loads(p.read_text()) if p.exists() else None


def _val(res, key):
    v = None if res is None else res.get(key)
    return v if v is not None and math.isfinite(v) else None


def cmd_R(args) -> int:
    """Determinism check only (amendment): |log rerun - log original| per run
    for dtmix_0.50/0.75/0.95 and whether it is exactly 0.  No H_rep verdict."""
    rows = []
    for r in R_RPMS:
        rerun_id = f"rec_l8_rpm{r}_th7"
        orig_id = next((c for c in (f"kmix_l8_rpm{r}", f"t5_l8_rpm{r}_th7")
                        if (ROOT / "runs" / c / "results.json").exists()), None)
        new, old = _res(rerun_id), (_res(orig_id) if orig_id else None)
        row = {"rpm": float(r), "rerun": rerun_id, "original": orig_id, "keys": {}}
        for k in DT_KEYS:
            d = _dlog(_val(old, k), _val(new, k))
            row["keys"][k] = {"rerun": _val(new, k), "original": _val(old, k),
                              "abs_log_diff": d,
                              "exactly_zero": None if d is None else bool(d == 0.0)}
        row["status"] = "ok" if all(v["abs_log_diff"] is not None
                                    for v in row["keys"].values()) else "missing"
        rows.append(row)
    n_ok = sum(r["status"] == "ok" for r in rows)
    n_zero = sum(all(v["exactly_zero"] for v in r["keys"].values())
                 for r in rows if r["status"] == "ok")
    res = {"n_reruns_available": n_ok, "n_bitwise_identical_dtmix": n_zero, "rows": rows}
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "test_r_rerun.json").write_text(json.dumps(_json_safe(res), indent=1))
    print(f"test R (determinism): {n_ok}/8 reruns available, "
          f"{n_zero} with all three dtmix exactly equal")
    return 0


# --------------------------------------------------------------- test P
P_BASES = {"22.5": ("kmix_l8_rpm22.5", "t5_l8_rpm22.4_th7", "t5_l8_rpm22.6_th7"),
           "30": ("kmix_l8_rpm30", "t5_l8_rpm29.9_th7", "t5_l8_rpm30.1_th7")}
P_OFFSETS = [("1e-2", 1e-2), ("1e-3", 1e-3), ("1e-6", 1e-6)]


def test_p_reading(med: dict) -> dict:
    """Pre-registered reading on medians (over bases) of |Delta log dtmix_0.95|
    per offset, keys 0.1, 1e-2, 1e-3, 1e-6 (None if not available).

    NOISE if median at 1e-6 >= 0.10; SMOOTH if at 1e-3 <= 0.02 AND the jumps
    fall by >= 10x per 2 decades.  Implementation of the second clause (spec
    silent): med(0.1)/med(1e-3) >= 10 (2 decades) and med(1e-2)/med(1e-6) >= 100
    (4 decades); a zero denominator counts as satisfied.  Otherwise mixed."""
    if any(med.get(k) is None for k in (0.1, 1e-2, 1e-3, 1e-6)):
        return {"verdict": "incomplete", "detail": "need all four offsets"}
    def fall(a, b, f):
        return b == 0 or a / b >= f
    if med[1e-6] >= 0.10:
        v = "NOISE"
    elif med[1e-3] <= 0.02 and fall(med[0.1], med[1e-3], 10) and fall(med[1e-2], med[1e-6], 100):
        v = "SMOOTH"
    else:
        v = "mixed"
    return {"verdict": v}


def plot_test_p(rows, path: Path) -> bool:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scripts import figstyle as fs
    plt.rcParams.update(fs.rcparams())
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    drew = False
    styles = {"22.5": ("#0072B2", "o"), "30": ("#D55E00", "s")}
    for base, (col, mk) in styles.items():
        pts = [(r["delta_rpm"], r["keys"]["dtmix_0.95"]["abs_log_diff"])
               for r in rows if r["base"] == base and r["keys"]["dtmix_0.95"]["abs_log_diff"] is not None]
        if pts:
            x, y = zip(*sorted(pts))
            y = [max(v, 1e-8) for v in y]      # exact zeros drawn at the floor
            ax.plot(x, y, color=col, marker=mk, ls=":", lw=1.1, label=f"{base} rpm")
            drew = True
    if drew:
        xs = np.array([1e-6, 1e-1])
        ax.plot(xs, 0.25 * xs / 0.1, color=fs.NEUTRAL, ls="--", lw=0.9, label="slope 1")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("rpm offset (rpm)")
        ax.set_ylabel("$|\\Delta \\ln\\, dt_{mix,0.95}|$ (-)")
        ax.grid(**fs.GRID_KW)
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
        fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return drew


def cmd_P(args) -> int:
    rows = []
    for base, (base_id, lo_id, hi_id) in P_BASES.items():
        b = _res(base_id)
        runs = [(f"p8_rpm{base}_d{tag}_th7", off) for tag, off in P_OFFSETS]
        runs.append((hi_id, 0.1))               # +0.1 side matches the sign of p8 offsets
        for rid, off in runs:
            r = _res(rid)
            row = {"base": base, "run": rid, "base_run": base_id, "delta_rpm": off, "keys": {}}
            for k in DT_KEYS + KLA_KEYS:
                row["keys"][k] = {"run": _val(r, k), "base": _val(b, k),
                                  "abs_log_diff": _dlog(_val(b, k), _val(r, k))}
            row["status"] = "ok" if row["keys"]["dtmix_0.95"]["abs_log_diff"] is not None else "missing"
            rows.append(row)
        r = _res(lo_id)                          # -0.1 side, reported (not used in the reading)
        row = {"base": base, "run": lo_id, "base_run": base_id, "delta_rpm": 0.1,
               "side": "-0.1", "keys": {}}
        for k in DT_KEYS + KLA_KEYS:
            row["keys"][k] = {"run": _val(r, k), "base": _val(b, k),
                              "abs_log_diff": _dlog(_val(b, k), _val(r, k))}
        row["status"] = "ok" if row["keys"]["dtmix_0.95"]["abs_log_diff"] is not None else "missing"
        rows.append(row)
    med = {}
    for off in (0.1, 1e-2, 1e-3, 1e-6):
        d = [r["keys"]["dtmix_0.95"]["abs_log_diff"] for r in rows
             if r["delta_rpm"] == off and r.get("side") != "-0.1"
             and r["keys"]["dtmix_0.95"]["abs_log_diff"] is not None]
        med[off] = float(np.median(d)) if len(d) == 2 else None
    reading = test_p_reading(med)
    res = {"pre_registered_reading": (
               "NOISE if the median |Delta| over the two bases at 1e-6 rpm >= 0.10; SMOOTH if at "
               "1e-3 rpm it is <= 0.02 AND the jumps fall by >= 10x per 2 decades; otherwise mixed."),
           "median_abs_dlog_dtmix_0.95_by_offset": {str(k): v for k, v in med.items()},
           "reading": reading, "rows": rows}
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "test_p_perturb.json").write_text(json.dumps(_json_safe(res), indent=1))
    plot_test_p(rows, OUT_DIR / "test_p_perturb.png")
    print("test P:", sum(r["status"] == "ok" for r in rows), "of", len(rows), "rows available;",
          "reading:", reading["verdict"])
    return 0


def cmd_check(args) -> int:
    s = chi_series(args.run_id)
    if s["status"] != "ok":
        print(args.run_id, s["status"], s.get("detail", ""))
        return 0
    print(f"{args.run_id}: {s['n_snapshots']} snapshots, window to {s['t_end_s']:.2f} s after release; "
          f"release offset {s['release_offset_s']:.3f} s (chi {s['release_offset_chi']:.4f})")
    print(json.dumps(_json_safe(compare_to_results(s)), indent=1))
    # In-progress runs have no results.json: compare against the solver's own chi(t)
    # rebuilt from tr_oxy.dat (same code path as postprocess), pointwise and by crossing.
    loc = locate(args.run_id)
    sol = solver_chi_from_dat(loc["data_dir"])
    if sol is not None:
        t, chi, k0, T_bio = sol
        t_s = (t - t[k0]) * T_bio
        ts = s["t_s"] + 0.0
        # snapshot times are relative to the same origin (t0 = first nonzero row)
        chi_sol = np.interp(ts, t_s[k0:], chi[k0:])
        m = (ts >= 0) & (ts <= t_s[-1])
        dchi = np.abs(chi_sol[m] - s["chi"]["native"][m])
        print(f"pointwise |chi_snap - chi_solver(tr_oxy)| over {int(m.sum())} snapshots: "
              f"max {dchi.max():.4g}, median {np.median(dchi):.4g}; solver chi range to "
              f"{chi[-1]:.3f} at {t_s[-1]:.2f} s, snapshots chi to {s['chi']['native'][m][-1]:.3f}")
        for thr in THRESHOLDS:
            a = crossing_time(t_s[k0:], chi[k0:], thr) if chi[-1] >= thr else math.nan
            b = s["dtmix"]["native"][f"{thr:.2f}"]
            print(f"  chi {thr:.2f}: solver tr_oxy interp {a:.4f} s, snapshots {b:.4f} s",
                  "" if not (math.isfinite(a) and math.isfinite(b)) else f"rel {b / a - 1:+.4%}")
    print("interface treatment:", compare_interface_treatment(args.run_id))
    print("dtmix by scale (s):")
    for k, v in s["dtmix"].items():
        print(f"  {k:8s}", {a: (None if not math.isfinite(b) else round(b, 3)) for a, b in v.items()})
    return 0



# ------------------------------------------------ solver dtmix quantisation
def solver_chi_from_dat(run_dir: Path):
    """The solver's own chi(t) rebuilt from tr_oxy.dat exactly as
    postprocess._compute_mixing_metrics does.  Returns (t_nd, chi, t0_index,
    T_bio) or None."""
    from scripts import postprocess as pp
    params = json.loads((run_dir / "params.json").read_text())
    try:
        t, c_sum, c_sum2 = pp._joined_cols(
            run_dir, "tr_oxy.dat", [pp._COL_C2_LIQ_SUM, pp._COL_C2_LIQ_SUM2])
        _, f_liq = pp._joined_cols(run_dir, "vol_frac_interf.dat", [pp._COL_F_LIQ])
    except (FileNotFoundError, ValueError):
        return None
    fm = float(f_liq.mean())
    sig = c_sum2 / fm - (c_sum / fm) ** 2
    nz = np.where(c_sum > 1e-10 * fm)[0]
    if not nz.size or sig[nz[0]] <= 0:
        return None
    k0 = int(nz[0])
    chi = np.clip(1.0 - sig / sig[k0], 0.0, 1.0)
    chi[:k0] = 0.0
    return t, chi, k0, pp._t_scales(params)[0]


def cmd_quant(args) -> int:
    """t_out cadence of tr_oxy.dat and the time quantisation of the solver's
    dtmix (first logged row with chi >= thr, no interpolation) against the same
    chi(t) linearly interpolated to the threshold, over kmix_l8_* and t5_l8_*."""
    rows = []
    for d in sorted(list((ROOT / "runs").glob("kmix_l8_*")) + list((ROOT / "runs").glob("t5_l8_*"))):
        r = solver_chi_from_dat(d)
        if r is None:
            continue
        t, chi, k0, T_bio = r
        dt_out = float(np.median(np.diff(t[k0:k0 + 200])))
        row = {"run": d.name, "rpm": json.loads((d / "params.json").read_text())["omega_b"] * 30 / math.pi,
               "t_out_nd": dt_out, "t_out_s": dt_out * T_bio, "dtmix": {}}
        for thr in THRESHOLDS:
            tt, cc = t[k0:] - t[k0], chi[k0:]
            first = crossing_time(tt, cc, thr, True) * T_bio
            interp = crossing_time(tt, cc, thr) * T_bio
            res = _res(d.name)
            row["dtmix"][f"{thr:.2f}"] = {
                "first_sample_s": first, "interp_s": interp,
                "results_json_s": _val(res, f"dtmix_{thr:.2f}"),
                "diff_s": first - interp,
                "diff_frac": (first - interp) / first if first and math.isfinite(first) else None}
        rows.append(row)
    summ = {}
    for thr in THRESHOLDS:
        k = f"{thr:.2f}"
        dd = [(abs(r["dtmix"][k]["diff_s"]), r["dtmix"][k]["diff_frac"], r["run"])
              for r in rows if r["dtmix"][k]["diff_frac"] is not None and math.isfinite(r["dtmix"][k]["diff_s"])]
        if dd:
            summ[k] = {"n_runs": len(dd),
                       "max_abs_diff_s": max(x[0] for x in dd),
                       "max_abs_diff_frac": max(abs(x[1]) for x in dd),
                       "median_abs_diff_frac": float(np.median([abs(x[1]) for x in dd])),
                       "run_of_max": max(dd, key=lambda x: x[0])[2]}
    res = {"t_out_s_by_run": {r["run"]: r["t_out_s"] for r in rows}, "summary": summ, "rows": rows}
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "solver_dtmix_quantisation.json").write_text(json.dumps(_json_safe(res), indent=1))
    print(json.dumps(_json_safe(summ), indent=1))
    ts = sorted({round(r["t_out_s"], 4) for r in rows})
    print("t_out_s values:", ts[:6], "...", ts[-3:], " t_out_nd values:", sorted({round(r['t_out_nd'], 4) for r in rows}))
    return 0


# --------------------------------------------------------------- test G
G_LIQ = 0.999         # liquid cell: f >= G_LIQ and cs >= G_LIQ
G_L_BIO_M = 0.25      # dx [m] = G_L_BIO_M / n  (coordinator: 0.25 m / 2^L)


def gradient_scale(f, cs, c) -> dict:
    """Scalar microscale lambda_c = sqrt(<c'^2> / <|grad c|^2>) over liquid
    cells, in cell units, and the fraction of liquid cells with 0.1 < c < 0.9.

    Mask (stated choice): a cell is liquid if f >= 0.999 and cs >= 0.999 (pure
    liquid, not an embedded cut cell).  Central differences
    (c[i+1]-c[i-1])/(2 dx), same in j, need the four neighbours: a cell enters
    the gradient average only if it and its four neighbours are liquid (and it
    is not on the domain edge).  c' = c - <c> with <c> the mean over the
    gradient-mask cells (unweighted; all cells are full cells).  The fraction
    0.1 < c < 0.9 uses the liquid cells (cell only, no stencil condition)."""
    liq = (f >= G_LIQ) & (cs >= G_LIQ)
    m = liq.copy()
    m[1:-1, 1:-1] &= liq[2:, 1:-1] & liq[:-2, 1:-1] & liq[1:-1, 2:] & liq[1:-1, :-2]
    m[0, :] = m[-1, :] = False
    m[:, 0] = m[:, -1] = False
    gx = np.zeros_like(c); gy = np.zeros_like(c)
    gx[1:-1, 1:-1] = (c[1:-1, 2:] - c[1:-1, :-2]) / 2.0
    gy[1:-1, 1:-1] = (c[2:, 1:-1] - c[:-2, 1:-1]) / 2.0
    g2 = (gx * gx + gy * gy)[m].mean()
    var = c[m].var()
    nl = int(liq.sum())
    return {"lambda_dx": float(math.sqrt(var / g2)) if g2 > 0 else math.nan,
            "frac_mid": float(((c > 0.1) & (c < 0.9))[liq].sum() / nl) if nl else math.nan,
            "n_grad_cells": int(m.sum()), "n_liquid": nl}


def gradient_series(run_id: str) -> dict:
    loc = locate(run_id)
    if loc["fields_dir"] is None or loc["data_dir"] is None:
        return {"run_id": run_id, "status": "missing"}
    ref = _solver_reference(loc)
    rows, started = [], False
    for snap in iter_snapshots(loc["fields_dir"], ["f", "cs", "c2"]):
        f, cs, c = snap["f"], snap["cs"], snap["c2"]
        if not started:
            if float((cs * f * c).sum()) <= 1e-10 * float((cs * f).sum()):
                continue
            started = True
        g = gradient_scale(f, cs, c)
        chi = None
        if ref["sigma0"]:
            chi = float(np.clip(1 - sigma2_native(f, cs, c, ref["f_mean"]) / ref["sigma0"], 0, 1))
        n = snap["n"]
        rows.append({"t_s": (snap["t"] - ref["t0"]) * ref["T_bio"], "chi_native": chi,
                     "n": n, "lambda_dx": g["lambda_dx"],
                     "lambda_m": g["lambda_dx"] * G_L_BIO_M / n, "frac_mid": g["frac_mid"]})
    if not rows:
        return {"run_id": run_id, "status": "missing", "detail": "no snapshots with tracer"}
    return {"run_id": run_id, "status": "ok", "rows": rows}


def _window_stats(series, lo=0.5, hi=0.95):
    """Statistics of lambda_c over snapshots with native chi in [lo, hi]."""
    r = [x for x in series["rows"] if x["chi_native"] is not None and lo <= x["chi_native"] <= hi
         and math.isfinite(x["lambda_dx"])]
    if not r:
        return None
    ld = np.array([x["lambda_dx"] for x in r])
    lm = np.array([x["lambda_m"] for x in r])
    return {"n_snapshots": len(r), "t_from_s": r[0]["t_s"], "t_to_s": r[-1]["t_s"],
            "lambda_dx_median": float(np.median(ld)), "lambda_m_median": float(np.median(lm)),
            "lambda_rel_scatter": float(np.std(lm) / np.mean(lm)),
            "frac_mid_median": float(np.median([x["frac_mid"] for x in r]))}


def test_g_verdict(ws8, ws9) -> dict:
    """Pre-registered rule for one rpm.  ratio = lambda_c(L8)/lambda_c(L9) of the
    window medians (metres).  Scatter check: combined relative scatter
    sqrt(s8^2 + s9^2) must be < 0.3 (the margin quoted in the spec)."""
    if ws8 is None or ws9 is None:
        return {"verdict": "missing"}
    ratio = ws8["lambda_m_median"] / ws9["lambda_m_median"]
    scat = math.hypot(ws8["lambda_rel_scatter"], ws9["lambda_rel_scatter"])
    l8, l9 = ws8["lambda_dx_median"], ws9["lambda_dx_median"]
    grid = (1 <= l8 <= 3) and (1 <= l9 <= 3) and (1.5 <= ratio <= 2.5)
    notgrid = (0.8 <= ratio <= 1.25) and l8 >= 5 and l9 >= 5
    return {"ratio_L8_over_L9": ratio, "lambda_dx_L8": l8, "lambda_dx_L9": l9,
            "combined_rel_scatter": scat, "scatter_below_0.3": bool(scat < 0.3),
            "verdict": "GRID-CONTROLLED" if grid else ("NOT grid-controlled" if notgrid else "neither")}


def plot_test_g(res_series, path: Path) -> bool:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scripts import figstyle as fs
    plt.rcParams.update(fs.rcparams())
    fig, axes = plt.subplots(2, 2, figsize=(10, 6.2), sharex="col")
    drew = False
    for j, rpm in enumerate((32.5, 37.5)):
        for lvl in (8, 9):
            s = res_series.get(f"scal_l{lvl}_rpm{rpm:g}_th7")
            if not s or s["status"] != "ok":
                continue
            t = [x["t_s"] for x in s["rows"]]
            col = fs.level_colour(lvl)
            axes[0, j].plot(t, [x["lambda_dx"] for x in s["rows"]], color=col, ls=":", lw=1.4, label=f"L{lvl}")
            axes[1, j].plot(t, [x["frac_mid"] for x in s["rows"]], color=col, ls=":", lw=1.4, label=f"L{lvl}")
            drew = True
        axes[0, j].set_title(f"{rpm:g} rpm, 7$^\\circ$")
        axes[1, j].set_xlabel("time after release (s)")
    axes[0, 0].set_ylabel("$\\lambda_c / \\Delta x$ (-)")
    axes[1, 0].set_ylabel("fraction of liquid cells, $0.1<c<0.9$ (-)")
    for ax in axes.ravel():
        ax.grid(**fs.GRID_KW)
    if drew:
        h, l = axes[0, 1].get_legend_handles_labels()
        axes[0, 1].legend(h, l, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
        fig.tight_layout()
        fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return drew


def cmd_G(args) -> int:
    ids = [f"scal_l{l}_rpm{r:g}_th7" for r in (32.5, 37.5) for l in (8, 9)] + \
          [f"rec_l8_rpm{r}_th7" for r in R_RPMS]
    ser = {i: gradient_series(i) for i in ids}
    res = {"runs": {}, "verdict_by_rpm": {}}
    for i, s in ser.items():
        res["runs"][i] = {"status": s["status"],
                          "window_chi_0.5_0.95": _window_stats(s) if s["status"] == "ok" else None,
                          "series": s.get("rows")}
    for r in (32.5, 37.5):
        res["verdict_by_rpm"][f"{r:g}"] = test_g_verdict(
            res["runs"][f"scal_l8_rpm{r:g}_th7"]["window_chi_0.5_0.95"],
            res["runs"][f"scal_l9_rpm{r:g}_th7"]["window_chi_0.5_0.95"])
    v = {x["verdict"] for x in res["verdict_by_rpm"].values()}
    res["verdict"] = ("missing" if "missing" in v else
                      v.pop() if len(v) == 1 else "neither/mixed across rpm")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "test_g_gradient_scale.json").write_text(json.dumps(_json_safe(res), indent=1))
    plot_test_g(ser, OUT_DIR / "test_g_gradient_scale.png")
    for i, s in ser.items():
        print(i, s["status"], (res["runs"][i]["window_chi_0.5_0.95"] or ""))
    print("test G verdict:", res["verdict"])
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("L").set_defaults(fn=cmd_L)
    sub.add_parser("R").set_defaults(fn=cmd_R)
    sub.add_parser("P").set_defaults(fn=cmd_P)
    sub.add_parser("G").set_defaults(fn=cmd_G)
    sub.add_parser("quant").set_defaults(fn=cmd_quant)
    pc = sub.add_parser("check")
    pc.add_argument("run_id")
    pc.set_defaults(fn=cmd_check)
    args = ap.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    raise SystemExit(main())
