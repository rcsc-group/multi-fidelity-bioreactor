"""Finish a sweep point that ran out of clock, without waiting to be told.

`submit_fig9.py` sizes each point's t_end as MARGIN x Kim's own dtmix_0.95,
on the premise that a coarser grid mixes faster than his. The first completed
point falsified the premise where it matters. At 37.5 rpm:

    chi = 0.50   ours  7.96 s   Kim  11.20 s   0.71x
    chi = 0.75   ours 17.02 s   Kim  16.80 s   1.01x
    chi = 0.95   ours    >44.7  Kim  33.61 s   >1.33x  (never reached)

Numerical diffusion buys the bulk homogenisation and nothing in the tail, so a
point can run to completion, pass every invariant, and still report NaN for
the one series Fig 9 is about. Extrapolating that point's own tail puts its
true dtmix_0.95 near 66 s -- about 2x Kim, not 1.3x -- which is a walltime no
single job can hold at the low-rpm end anyway.

Rather than re-guess the margin, a point that stops short is continued from
its own checkpoint until it reaches the threshold. That is strictly cheaper
than a longer cold start (the 80-cycle spin-up is paid once), it keeps every
job inside the 48 h cap, and postprocess already stitches the segments
(`_parent_run`, see `_segment_chain`) so the joined chi starts at injection
however many segments it spans.

The one thing this must NOT do is spend a segment on a point whose NaN is not
a clock problem. A sigma^2_max away from its analytic 0.25 means chi is
normalised against the wrong reference; the threshold is then unreachable by
construction and more compute returns the same NaN a day later.

Usage:
    uv run python scripts/autoextend.py                # report, submit nothing
    uv run python scripts/autoextend.py --submit       # continue what is short
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
from scripts.postprocess import (           # noqa: E402
    _COL_C2_LIQ_SUM, _COL_C2_LIQ_SUM2, _COL_F_LIQ,
    _SIGMA2_MAX_EXACT, _SIGMA2_MAX_TOL, _joined_cols, _t_scales,
)

TARGET_KEY = "dtmix_0.95"
TARGET_CHI = 0.95
DEFAULT_MARGIN = 1.25    # the tail fit is an extrapolation, not a measurement
DEFAULT_TAIL_FRAC = 0.35
DEFAULT_MAX_CYCLES = 400.0
MAX_SEGMENTS = 3         # a point needing a 4th segment is a finding, not a job


# ── which points deserve more compute ────────────────────────────────────────

def needs_extension(results: dict, target_key: str = TARGET_KEY) -> bool:
    """True only when the threshold is unreached AND the normalisation is sound.

    A missing key means postprocess never reached the mixing block at all
    (no tr_oxy.dat, too few rows) -- that is a broken run, not a short one.
    """
    if target_key not in results:
        return False
    value = results[target_key]
    if value is not None and not (isinstance(value, float) and math.isnan(value)):
        return False
    sigma2 = results.get("sigma2_max")
    if sigma2 is None or (isinstance(sigma2, float) and math.isnan(sigma2)):
        return False
    lo = _SIGMA2_MAX_EXACT * (1 - _SIGMA2_MAX_TOL)
    hi = _SIGMA2_MAX_EXACT * (1 + _SIGMA2_MAX_TOL)
    return lo <= float(sigma2) <= hi


# ── how much more compute ────────────────────────────────────────────────────

def extension_cycles(t_s, chi, T_per: float, target: float = TARGET_CHI,
                     margin: float = DEFAULT_MARGIN,
                     tail_frac: float = DEFAULT_TAIL_FRAC,
                     max_cycles: float = DEFAULT_MAX_CYCLES) -> float:
    """Rocking cycles still needed to reach `target`, from the curve's own tail.

    The unmixed fraction u = 1 - chi decays exponentially with a time constant
    that is NOT constant: the measured tail relaxes several times slower than
    the bulk. Fitting the whole curve therefore under-estimates the remainder
    and re-submits a segment that falls short again, so the fit uses only the
    last `tail_frac` of the samples.

    `t_s` is in seconds and so is the fit; the return is in cycles because that
    is what a walltime is sized from.
    """
    t_s = np.asarray(t_s, dtype=float)
    chi = np.asarray(chi, dtype=float)
    u_end = 1.0 - float(chi[-1])
    u_target = 1.0 - target
    if u_end <= u_target:
        return 0.0

    n_tail = max(3, int(round(tail_frac * t_s.size)))
    t_fit, u_fit = t_s[-n_tail:], 1.0 - chi[-n_tail:]
    if np.any(u_fit <= 0) or np.ptp(t_fit) <= 0:
        return max_cycles
    slope = np.polyfit(t_fit, np.log(u_fit), 1)[0]
    if not np.isfinite(slope) or slope >= 0:
        return max_cycles          # flat or rising: no decay to extrapolate

    tau = -1.0 / slope
    extra_s = tau * math.log(u_end / u_target) * margin
    if not np.isfinite(extra_s) or extra_s <= 0:
        return max_cycles
    return float(min(extra_s / T_per, max_cycles))


def walltime_safety(confidence: str) -> float:
    """Cushion on the cost-model estimate, by how much that estimate is trusted.

    A "measured" cost is for this exact rpm and needs only the usual margin
    (the model is measured pre-injection; tracer gradients cost more). A
    fallback cost is the MAX over the rpms measured at that (level, ntasks) --
    which, with a single entry in the table, is no maximum at all. Cost per
    cycle rises steeply as rpm drops (2.9x from 32.5 to 17.5 rpm at L10/32),
    so a fallback must cover that spread or the segment times out and is lost.
    """
    return 1.6 if confidence == "measured" else 3.2


def period_seconds(params: dict) -> float:
    """One rocking period in SECONDS.

    Not `_t_scales(params)[1]` -- that is the period in the solver's
    non-dimensional units, and using it to convert an extrapolated number of
    seconds into cycles inflates the count by T_bio (2.63x at 37.5 rpm).
    """
    return 2.0 * math.pi / float(params["omega_b"])


def chi_from_run(run_dir: Path, params: dict):
    """(t in seconds since injection, chi) for the whole segment chain."""
    t_nd, c_sum, c_sum2 = _joined_cols(
        run_dir, "tr_oxy.dat", [_COL_C2_LIQ_SUM, _COL_C2_LIQ_SUM2])
    _, f_liq = _joined_cols(run_dir, "vol_frac_interf.dat", [_COL_F_LIQ])
    f_mean = float(f_liq.mean())
    c_mean, c_mean2 = c_sum / f_mean, c_sum2 / f_mean
    sigma2 = c_mean2 - c_mean ** 2
    nonzero = np.where(c_sum > 1e-10 * f_mean)[0]
    if nonzero.size == 0:
        raise ValueError(f"{run_dir.name}: tracer never injected")
    i0 = int(nonzero[0])
    sigma2_max = float(sigma2[i0])
    if sigma2_max <= 0:
        raise ValueError(f"{run_dir.name}: sigma2_max={sigma2_max}")
    chi = np.clip(1.0 - sigma2 / sigma2_max, 0.0, 1.0)[i0:]
    T_bio, _ = _t_scales(params)
    return (t_nd[i0:] - t_nd[i0]) * T_bio, chi


# ── where the next segment attaches ──────────────────────────────────────────

def chain_tip(runs_root: Path, base_run_id: str) -> str:
    """The deepest `<base>_ext<N>` that exists, else the base run itself.

    Matched with an anchored pattern: `fig9_l9_rpm3` is a plain-string prefix
    of `fig9_l9_rpm30`, and a prefix match would attach rpm30's continuation
    to rpm3's chain.
    """
    pattern = re.compile(rf"^{re.escape(base_run_id)}_ext(\d+)$")
    best, best_n = base_run_id, 0
    for d in Path(runs_root).iterdir():
        if not d.is_dir():
            continue
        m = pattern.match(d.name)
        if m and int(m.group(1)) > best_n:
            best, best_n = d.name, int(m.group(1))
    return best


def segment_index(run_id: str) -> int:
    """0 for a base run, N for `<base>_ext<N>`."""
    m = re.search(r"_ext(\d+)$", run_id)
    return int(m.group(1)) if m else 0


def tip_results(runs_root: Path, base_run_id: str) -> dict:
    """The finished measurement for a point: the LAST segment's results.json.

    A base run that was extended keeps its own NaN forever; anything reading
    the base run_id directly (plot_fig9) would report the point as missing.
    """
    runs_root = Path(runs_root)
    for name in (chain_tip(runs_root, base_run_id), base_run_id):
        rj = runs_root / name / "results.json"
        if rj.exists():
            return json.loads(rj.read_text())
    return {}


# ── submission ───────────────────────────────────────────────────────────────

def plan_extension(runs_root: Path, base_run_id: str,
                   margin: float = DEFAULT_MARGIN,
                   target_key: str = TARGET_KEY,
                   target_chi: float = TARGET_CHI) -> dict | None:
    """What the next segment for this point should be, or None if it needs none."""
    runs_root = Path(runs_root)
    tip = chain_tip(runs_root, base_run_id)
    tip_dir = runs_root / tip
    rj = tip_dir / "results.json"
    if not rj.exists():
        return None                       # still running; nothing to judge yet
    results = json.loads(rj.read_text())
    if not needs_extension(results, target_key=target_key):
        return None
    if segment_index(tip) >= MAX_SEGMENTS:
        return {"run_id": tip, "blocked": True, "reason":
                f"{MAX_SEGMENTS} segments already spent and chi is still short"}
    params = json.loads((tip_dir / "params.json").read_text())
    T_per = period_seconds(params)
    t_s, chi = chi_from_run(tip_dir, params)
    cycles = extension_cycles(t_s, chi, T_per=T_per, target=target_chi,
                              margin=margin)
    if cycles <= 0:
        return None
    return {"parent": tip, "run_id": f"{base_run_id}_ext{segment_index(tip) + 1}",
            "cycles": cycles, "T_per": T_per, "params": params,
            "chi_now": float(chi[-1]), "target_chi": target_chi,
            "blocked": False}


def submit_extension(plan: dict, ntasks: int, min_per_cycle: float,
                     walltime_safety: float = 1.6, dry: bool = False) -> str | None:
    from scripts.simulate import submit_slurm
    from scripts.dump_fields import fields as dump_fields

    runs_root, parent = ROOT / "runs", plan["parent"]
    ckpt = runs_root / parent / "checkpoint.dump"
    if not ckpt.exists():
        print(f"  {plan['run_id']}: no checkpoint in {parent}; cannot continue")
        return None
    # The restart branch is gated on params.t_checkpoint > 0, not on a dump
    # being staged -- without it the binary silently COLD starts (see
    # submit_fig9_l10_anchor.py).
    t_ck = dump_fields(ckpt)[0]["t"]
    if t_ck <= 0:
        print(f"  {plan['run_id']}: checkpoint reports t={t_ck}; cannot arm a restart")
        return None

    T_bio, T_per_nd = _t_scales(plan["params"])
    params = dict(plan["params"])
    params.update({
        "run_id": plan["run_id"],
        "t_checkpoint": t_ck,
        "t_end": round(plan["cycles"] * T_per_nd, 4),
        "n_mix_cycles": 1,          # continuation: the tracer is already in
        "restart_continue": 1,      # keep tracers, do not re-inject
        "_parent_run": parent,
    })
    hours = plan["cycles"] * min_per_cycle / 60.0 * walltime_safety
    walltime = f"{min(int(hours) + 1, 47):02d}:00:00"
    print(f"  {plan['run_id']}: chi={plan['chi_now']:.4f} -> "
          f"{plan.get('target_chi', TARGET_CHI)}, "
          f"{plan['cycles']:.0f} cyc, t_end={params['t_end']}, walltime={walltime}")
    if dry:
        print("       DRY RUN -- not submitted")
        return None
    job = submit_slurm(params, project_root=ROOT, runs_root=runs_root,
                       walltime=walltime,
                       template=ROOT / "config" / "slurm_mpi_template.sh",
                       cpus=1, ntasks=ntasks, mem="4G", checkpoint=str(ckpt))
    rec = ROOT / "logs" / f"submitted_{job}.txt"
    rec.parent.mkdir(exist_ok=True)
    rec.write_text(f"{plan['run_id']}\n")
    print(f"       submitted {plan['run_id']}  job={job}")
    return job


def main() -> None:
    from scripts.cost_model import min_per_cycle

    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", default="fig9_l9_rpm")
    ap.add_argument("--level", type=int, default=9)
    ap.add_argument("--ntasks", type=int, default=32)
    ap.add_argument("--submit", action="store_true",
                    help="without this the planner only reports")
    a = ap.parse_args()

    runs_root = ROOT / "runs"
    bases = sorted({d.name for d in runs_root.iterdir()
                    if d.is_dir() and d.name.startswith(a.prefix)
                    and "_ext" not in d.name})
    for base in bases:
        plan = plan_extension(runs_root, base)
        if plan is None:
            continue
        if plan.get("blocked"):
            print(f"  {base}: BLOCKED -- {plan['reason']}")
            continue
        mpc, conf = min_per_cycle(level=a.level, ntasks=a.ntasks,
                                  rpm=plan["params"]["omega_b"] * 60 / (2 * math.pi))
        submit_extension(plan, ntasks=a.ntasks, min_per_cycle=mpc,
                         walltime_safety=walltime_safety(conf),
                         dry=not a.submit)


if __name__ == "__main__":
    main()
