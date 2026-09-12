"""Settling-time model: settle_cycle(run, ref, rpm_to, theta_to, level) ->
first cycle (from the run's own checkpoint time) after which the run's
tau_mean waveform stays within tol of an INDEPENDENT reference's converged
waveform, for the rest of the run.

This generalizes analyze_settling_aligned.py's validated metric (RMS
distance between each cycle's phase-resampled waveform and a converged
target waveform, with the run's waveform cross-correlation-aligned against
the target to remove analysis-side phase error -- see diary.md 2026-09-08)
to an arbitrary (run, reference) pair and (rpm, theta) condition, instead of
one hardcoded run at one hardcoded condition.

Deliberately NOT the per-cycle-peak metric used in settling_study_v2 /
analyze_settling_after_fix.py -- that metric was proven void on 2026-09-08
(it tracks the 3-cycle forcing ramp, not physical convergence: alpha(t)
reaches 1.0 at the end of cycle 2 regardless of the true transient length,
so a per-cycle MAX always reads "converged" by cycle 2).

Usage as a library:
    from scripts.settling_model import settle_cycle
    k = settle_cycle(run_dir, ref_dir, rpm_to=32.5, theta_to=7.0, tol=0.02)

Usage as a script (recompute the settling grid from existing raw data,
zero new compute -- all 24 transition/reference run pairs from the
settling_study_v2 manifest are still on disk):
    uv run python scripts/settling_model.py
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent.parent
RUNS = ROOT / "runs"
NPH = 256
_PHASE = np.linspace(0, 1, NPH, endpoint=False)


def t_per_nd(rpm: float, theta: float) -> float:
    """Non-dimensional rocking period T_per / T_bio for this (rpm, theta)."""
    w = rpm * 2 * math.pi / 60.0
    L, H = 0.25, 2 * 0.03575
    T = 2 * math.pi / w
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(theta)))
    U = V / (H * 0.5) / T
    T_bio = L / U
    return T / T_bio


def waveforms_from_series(t: np.ndarray, y: np.ndarray, T_nd: float, t0: float = 0.0, nph: int = NPH) -> list[np.ndarray | None]:
    """Per-cycle waveforms from raw (t, y) arrays, phase-resampled onto a
    common [0,1) grid. Cycle k covers non-dimensional time
    [t0 + k*T_nd, t0 + (k+1)*T_nd). A cycle with fewer than 8 samples (can't
    resolve a waveform shape) is returned as None rather than silently
    interpolated from noise. Used directly by waveforms() (single run) and
    by analyze_extended_settling.py (a run concatenated with its own later
    same-condition extension, to see past what one run's own recorded
    window could show)."""
    c = (t - t0) / T_nd
    n_cycles = int(c[-1]) if len(c) else 0
    phase = np.linspace(0, 1, nph, endpoint=False)
    out: list[np.ndarray | None] = []
    for k in range(n_cycles):
        m = (c >= k) & (c < k + 1)
        if m.sum() < 8:
            out.append(None)
            continue
        ph, v = c[m] - k, y[m]
        out.append(np.interp(phase, ph, v, period=1.0))
    return out


def waveforms(run_dir: Path, T_nd: float, t0: float = 0.0, nph: int = NPH) -> list[np.ndarray | None]:
    """Per-cycle waveforms for one run's shear_stress.dat -- see
    waveforms_from_series() for the resampling logic."""
    d = np.loadtxt(run_dir / "shear_stress.dat", skiprows=1)
    return waveforms_from_series(d[:, 1], d[:, 5], T_nd, t0=t0, nph=nph)


def converged_waveform(waves: list[np.ndarray | None], n_tail: int = 20) -> np.ndarray | None:
    """Mean waveform over the last n_tail available cycles (None entries skipped)."""
    tail = [w for w in waves[-n_tail:] if w is not None]
    if not tail:
        return None
    return np.mean(tail, axis=0)


def best_phase_shift(cand_conv: np.ndarray, target_conv: np.ndarray, nph: int = NPH) -> float:
    """Sub-sample phase shift (in cycle fraction, [-0.5, 0.5)) that minimizes
    the RMS residual between cand_conv and target_conv. Removes constant
    analysis-side phase error (see diary.md 2026-09-08, trap #1) before any
    residual is trusted as a physical transient."""
    phase = np.linspace(0, 1, nph, endpoint=False)
    grid = np.linspace(-0.5, 0.5, 2001)
    best, bshift = None, 0.0
    for s in grid:
        q = np.mod(phase + s, 1.0)
        cand = np.interp(q, phase, cand_conv, period=1.0)
        r = np.sqrt(np.mean((cand - target_conv) ** 2))
        if best is None or r < best:
            best, bshift = r, s
    return bshift


def residual_series(
    t: np.ndarray, y: np.ndarray, ref_dir: Path, rpm_to: float, theta_to: float,
    t_checkpoint: float = 0.0, n_tail: int = 20, nph: int = NPH,
) -> tuple[np.ndarray | None, float]:
    """Per-cycle RMS residual of (t, y) against ref_dir's own converged
    waveform, phase-aligned (see module docstring). Returns (resid, scale);
    resid is None if (t, y) has no usable converged tail of its own to align
    from. scale is RMS(ref's converged waveform) -- the normalization every
    residual and every noise-floor number in this module is expressed
    relative to.
    """
    T = t_per_nd(rpm_to, theta_to)
    ref_waves = waveforms(ref_dir, T, t0=0.0, nph=nph)
    ref_conv = converged_waveform(ref_waves, n_tail)
    if ref_conv is None:
        raise ValueError(f"reference {ref_dir} has no usable converged waveform")
    scale = float(np.sqrt(np.mean(ref_conv ** 2)))

    run_waves = waveforms_from_series(t, y, T, t0=t_checkpoint, nph=nph)
    run_conv = converged_waveform(run_waves, n_tail)
    if run_conv is None:
        return None, scale
    shift = best_phase_shift(run_conv, ref_conv, nph)

    phase = np.linspace(0, 1, nph, endpoint=False)
    q = np.mod(phase + shift, 1.0)
    resid = np.full(len(run_waves), np.nan)
    for k, w in enumerate(run_waves):
        if w is None:
            continue
        w_aligned = np.interp(q, phase, w, period=1.0)
        resid[k] = np.sqrt(np.mean((w_aligned - ref_conv) ** 2)) / scale
    return resid, scale


def _settle_from_residual(resid: np.ndarray | None, tol: float) -> int | None:
    if resid is None:
        return None
    for k in range(len(resid)):
        seg = resid[k:]
        valid = seg[~np.isnan(seg)]
        if len(valid) > 0 and np.all(valid < tol):
            return k
    return None


def settle_cycle(
    run_dir: Path,
    ref_dir: Path,
    rpm_to: float,
    theta_to: float,
    t_checkpoint: float = 0.0,
    tol: float = 0.02,
    n_tail: int = 20,
    nph: int = NPH,
) -> tuple[int | None, int]:
    """First cycle k (0-based, from t_checkpoint) after which run_dir's
    waveform stays within tol*scale of ref_dir's own converged waveform for
    every remaining cycle. Returns (k, n_cycles_available); k is None if it
    never settles within the run's available data (NOT the same as "settled
    at the last cycle" -- callers must not confuse the two).

    ref_dir must be an INDEPENDENT run at (rpm_to, theta_to) that has itself
    reached a converged quasi-steady state (its own last n_tail cycles are
    stable) -- this is the ground truth the run is judged against, not a
    self-referential check.

    NOTE: a fixed tol is a known limitation -- see reference_noise_floor()
    and settle_cycle_adaptive() below, added after 2026-09-12's finding that
    different (rpm, theta) conditions have genuinely different natural
    fluctuation floors (0.8% to 5.5% observed), so a flat 2% both
    false-negatives calm conditions that have actually converged and would
    false-positive if applied to a noisier condition mid-transient.
    """
    d = np.loadtxt(run_dir / "shear_stress.dat", skiprows=1)
    resid, _ = residual_series(d[:, 1], d[:, 5], ref_dir, rpm_to, theta_to,
                                t_checkpoint=t_checkpoint, n_tail=n_tail, nph=nph)
    n_avail = int((d[-1, 1] - t_checkpoint) / t_per_nd(rpm_to, theta_to)) if len(d) else 0
    return _settle_from_residual(resid, tol), n_avail


def reference_noise_floor(
    ref_dir: Path, rpm: float, theta: float, n_tail: int = 20, nph: int = NPH,
) -> tuple[float, float]:
    """(mean, max) of a reference run's own cycle-to-cycle RMS residual
    against ITS OWN converged (last-n_tail-cycle) waveform -- i.e. how much
    this condition's quasi-steady state naturally fluctuates on its own,
    independent of any transient. Found 2026-09-12 to range 0.8%-5.5% across
    the settling-study grid: a single fixed settle tolerance cannot be
    correct for all conditions, so this is the per-condition scale to
    compare a transient's residual against instead.
    """
    T = t_per_nd(rpm, theta)
    waves = waveforms(ref_dir, T, nph=nph)
    tail = [w for w in waves[-n_tail:] if w is not None]
    if not tail:
        raise ValueError(f"reference {ref_dir} has no usable tail")
    conv = np.mean(tail, axis=0)
    scale = np.sqrt(np.mean(conv ** 2))
    phase = np.linspace(0, 1, nph, endpoint=False)
    resids = []
    for w in tail:
        shift = best_phase_shift(w, conv, nph)
        q = np.mod(phase + shift, 1.0)
        w_aligned = np.interp(q, phase, w, period=1.0)
        resids.append(np.sqrt(np.mean((w_aligned - conv) ** 2)) / scale)
    return float(np.mean(resids)), float(np.max(resids))


def adaptive_tol(ref_dir: Path, rpm: float, theta: float, margin: float = 1.5,
                  min_tol: float = 0.02, n_tail: int = 20) -> float:
    """Per-condition settle tolerance: margin x this condition's own natural
    noise ceiling (reference_noise_floor's max, not mean -- a transient
    should be judged converged only once it's within the reference's own
    worst cycle-to-cycle wobble, not just its average), floored at min_tol so
    an unusually quiet condition doesn't end up chasing sub-percent noise.
    """
    _, ref_max = reference_noise_floor(ref_dir, rpm, theta, n_tail=n_tail)
    return max(min_tol, margin * ref_max)


def settle_cycle_adaptive(
    run_dir: Path, ref_dir: Path, rpm_to: float, theta_to: float,
    t_checkpoint: float = 0.0, margin: float = 1.5, min_tol: float = 0.02,
    n_tail: int = 20, nph: int = NPH,
) -> tuple[int | None, int, float]:
    """Like settle_cycle(), but tol is derived from the reference's own
    noise floor instead of a fixed number. Returns (k, n_cycles_available,
    tol_used)."""
    tol = adaptive_tol(ref_dir, rpm_to, theta_to, margin=margin, min_tol=min_tol, n_tail=n_tail)
    k, n_avail = settle_cycle(run_dir, ref_dir, rpm_to, theta_to,
                               t_checkpoint=t_checkpoint, tol=tol, n_tail=n_tail, nph=nph)
    return k, n_avail, tol


# NOTE on a self-referential ("no independent reference") drift check:
# tried on 2026-09-12 (comparing waveform SHAPE, not just per-cycle peak,
# across tail sub-windows, hoping to catch what the peak-based check in
# validate_run.py's check_drift() misses). Tested directly against
# 7bb073d2 (KNOWN, via an independent reference, to still be ~15-20% from
# converged at 60 cycles) and 7e103866 (KNOWN converged): both gave a
# similar drift/local-noise ratio (~0.8-0.9). A slowly-decaying signal's
# tail is, by construction, locally indistinguishable from noise around a
# not-yet-converged mean -- no purely self-referential tail statistic can
# reliably tell the two cases apart. Do not add one here expecting it to
# work; the real check requires an independent reference
# (settle_cycle_adaptive, above). Where no reference exists, the honest
# defense is a well-margined a-priori settling_cycles() estimate plus
# extension, not an after-the-fact certifier.


_FIT_CACHE: list[dict] | None = None


def _load_fit() -> list[dict]:
    global _FIT_CACHE
    if _FIT_CACHE is None:
        path = ROOT / "experiments" / "settling_model_fit.json"
        _FIT_CACHE = json.loads(path.read_text()) if path.exists() else []
    return _FIT_CACHE


class UnresolvedSettlingCondition(Exception):
    """Raised when the requested transition matches a condition known to be
    a flagged anomaly (never settles to the reference within any tested
    window) -- never silently return a number for this, it would be a guess
    dressed up as data. Caller must get a human decision."""


def settling_cycles(level: int, from_rpm: float, from_theta: float,
                     to_rpm: float, to_theta: float,
                     safety_margin: float = 1.5) -> tuple[int, str]:
    """Planned transition cycle count for chain.py's auto-sizing, with a
    safety margin applied. Returns (cycles, confidence):
      "measured"           -- exact (level, from, to) match with a real
                               settle_cycle from fit_settling_model.py.
      "low_confidence_max"  -- no exact match (new level/condition, or a
                               specific edge that stayed censored even after
                               extension); falls back to the largest
                               MEASURED settle_cycle anywhere in the fitted
                               table x an extra margin, since we've directly
                               observed settling can take 100+ cycles and a
                               small default would repeat this project's
                               core mistake.

    Raises UnresolvedSettlingCondition for a transition matching a known
    flagged anomaly (2026-09-12: 32.5rpm theta 7->2, all levels) -- that
    condition does not have a settling time to estimate; it needs a human
    decision, not an extrapolated number.
    """
    rows = _load_fit()
    if not rows:
        raise RuntimeError("settling_model_fit.json not found -- run "
                            "scripts/fit_settling_model.py first")

    for r in rows:
        if r["confidence"] == "flagged_anomaly" and \
           abs(r["from"]["rpm"] - from_rpm) < 1e-6 and abs(r["from"]["theta"] - from_theta) < 1e-6 and \
           abs(r["to"]["rpm"] - to_rpm) < 1e-6 and abs(r["to"]["theta"] - to_theta) < 1e-6:
            raise UnresolvedSettlingCondition(
                f"L{level} {from_rpm:g}/{from_theta:g} -> {to_rpm:g}/{to_theta:g} is a "
                f"flagged anomaly (see {r.get('note')}) -- do not auto-size this transition")

    for r in rows:
        if r["confidence"] == "measured" and r["level"] == level and \
           abs(r["from"]["rpm"] - from_rpm) < 1e-6 and abs(r["from"]["theta"] - from_theta) < 1e-6 and \
           abs(r["to"]["rpm"] - to_rpm) < 1e-6 and abs(r["to"]["theta"] - to_theta) < 1e-6:
            return math.ceil(r["settle_cycle"] * safety_margin), "measured"

    measured = [r["settle_cycle"] for r in rows if r["confidence"] == "measured"]
    if not measured:
        raise RuntimeError("no measured settling data at all -- cannot even fall back conservatively")
    fallback = math.ceil(max(measured) * safety_margin * 1.3)  # extra margin on top: this is an
    # UNTESTED condition, not just an untested exact edge with siblings nearby
    return fallback, "low_confidence_max"


def _recompute_settling_grid(tol: float = 0.02) -> list[dict]:
    """Recompute settle_cycle for every edge in the settling_study_v2
    manifest, using raw shear_stress.dat still on disk -- zero new compute.
    """
    manifest = json.loads((ROOT / "experiments" / "settling_study_v2_manifest.json").read_text())
    results = []
    for level_str, edges in manifest["edges"].items():
        level = int(level_str)
        for tag, e in edges.items():
            frm, to, run_id = e["from"], e["to"], e["run_id"]
            ref_id = f"settling_ref_L{level}_rpm{to['rpm']:g}_th{to['theta']:g}"
            params = json.loads((RUNS / run_id / "params.json").read_text())
            t_ckpt = params["t_checkpoint"]
            k, n_avail = settle_cycle(
                RUNS / run_id, RUNS / ref_id,
                rpm_to=to["rpm"], theta_to=to["theta"],
                t_checkpoint=t_ckpt, tol=tol,
            )
            results.append({
                "level": level, "from": frm, "to": to,
                "d_rpm": to["rpm"] - frm["rpm"], "d_theta": to["theta"] - frm["theta"],
                "settle_cycle": k, "n_cycles_available": n_avail, "tol": tol,
            })
    return results


if __name__ == "__main__":
    for tol in (0.05, 0.02, 0.01):
        rows = _recompute_settling_grid(tol=tol)
        never = sum(1 for r in rows if r["settle_cycle"] is None)
        print(f"\n=== tol={tol:.0%} ({never}/{len(rows)} never settled within available data) ===")
        print(f"{'L':>2} {'from':>14} {'to':>14} {'settle_cycle':>12} {'n_avail':>8}")
        for r in rows:
            frm_s = f"{r['from']['rpm']:g}/{r['from']['theta']:g}"
            to_s = f"{r['to']['rpm']:g}/{r['to']['theta']:g}"
            sc = "NEVER" if r["settle_cycle"] is None else r["settle_cycle"]
            print(f"{r['level']:>2} {frm_s:>14} {to_s:>14} {sc!s:>12} {r['n_cycles_available']:>8}")
        out = ROOT / "experiments" / f"settling_model_grid_tol{int(tol*100)}.json"
        out.write_text(json.dumps(rows, indent=2))
        print(f"saved {out}")
