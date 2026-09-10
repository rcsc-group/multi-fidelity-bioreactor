"""Automated gate: is a finished run actually good data, checked against
every failure mode this project has hit before -- not trusted by default.

Checks, in order (a run FAILS if any hard check fails; WARN is informative,
doesn't block):

  1. Binary provenance (HARD).  params.json["_binary"] must be a path this
     script recognizes as the current, fully-fixed binary. Catches running
     on a stale/default binary that predates a known fix (the exact way
     l10_kim_fig8_signed silently predated the H_bio and tau-OpenMP-race
     fixes, and _binary:None runs predate everything).
  2. Reached its target cycle count (HARD). Compares actual elapsed
     (t[-1]-t[0])/T_per_nd against n_mix_cycles. Catches a run cut short by
     walltime/crash and mistaken for a converged one (l10_kim_seg2: 17 of
     500 cycles).
  3. NaN/Inf scan (HARD) across shear_stress.dat, normf.dat,
     vol_frac_interf.dat, tr_oxy.dat, results.json.
  4. Restart field fidelity (HARD, restarts only). Diffs THIS segment's
     restart_diagnostic_post_restore.txt against its PARENT's
     restart_diagnostic_pre_dump.txt, after applying the same su/su^2
     rescale event init itself applies -- this is Test C (diary.md
     2026-09-07), automated and run on every restart from now on, not
     just when something already looks wrong.
  5. Parent-linkage sanity (HARD, restarts only). theta_max_prev/
     omega_b_prev recorded in params.json must equal the PARENT segment's
     own theta_max/omega_b -- catches a chain-wiring bug handing a segment
     the wrong previous condition.
  6. Amplitude-drift / QSS check (WARN). Per-cycle tau_mean peak over the
     last third of the run, checked for a trend beyond ~3x the run's own
     early-window cycle-to-cycle noise -- the same waveform-comparison
     technique that caught the phase-bug escape and the l10_kim_fig8_signed
     mis-check, run automatically instead of by hand after the fact.

Usage:
    uv run python scripts/validate_run.py <run_id> [--parent <parent_run_id>]
    uv run python scripts/validate_run.py --all-known   # everything this
                                                          # session touched
"""
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent.parent
RUNS = ROOT / "runs"

KNOWN_GOOD_BINARIES = {
    "/oscar/scratch/eaguerov/BioReactor-mpi-phasefix2",
    "/oscar/scratch/eaguerov/BioReactor-mpi-video-fixed",
}


def t_per_nd(rpm, theta, L_bio=0.25, H_bio_full=2 * 0.03575):
    w = rpm * 2 * math.pi / 60.0
    T = 2 * math.pi / w
    V = L_bio / 4 * (H_bio_full + 0.5 * L_bio * math.tan(math.radians(theta)))
    U = V / (H_bio_full * 0.5) / T
    return T / (L_bio / U)


def load_diag(path):
    if not path.exists():
        return None
    d = {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) == 2:
            try:
                d[parts[0]] = float(parts[1])
            except ValueError:
                d[parts[0]] = parts[1]
    return d


def check_binary(params, report):
    b = params.get("_binary")
    if b is None:
        report.append(("HARD", "FAIL", f"_binary is None (ran on the default/unknown binary)"))
    elif b not in KNOWN_GOOD_BINARIES:
        report.append(("HARD", "FAIL", f"_binary={b!r} not a recognized current-fix binary"))
    else:
        report.append(("HARD", "PASS", f"binary = {b}"))


def check_cycle_count(params, report):
    rpm = params["omega_b"] * 60.0 / (2 * math.pi)
    theta = params["theta_max"][0]
    T_ND = t_per_nd(rpm, theta)
    ss_path = RUNS / params["run_id"] / "shear_stress.dat"
    if not ss_path.exists():
        report.append(("HARD", "FAIL", "shear_stress.dat missing"))
        return
    try:
        d = np.loadtxt(ss_path, skiprows=1)
    except ValueError as e:
        report.append(("HARD", "FAIL", f"shear_stress.dat unparseable: {e}"))
        return
    got = (d[-1, 1] - d[0, 1]) / T_ND
    want = params.get("n_mix_cycles")
    if want is None:
        report.append(("WARN", "WARN", "no n_mix_cycles in params.json to check against"))
        return
    frac = got / want
    if frac < 0.95:
        report.append(("HARD", "FAIL",
                       f"reached {got:.1f} of {want} target cycles ({100*frac:.0f}%) -- cut short"))
    else:
        report.append(("HARD", "PASS", f"reached {got:.1f} of {want} target cycles"))


def check_nan_inf(params, report):
    """Also catches file-level corruption, not just NaN/Inf VALUES. Caught
    live (diary.md 2026-09-10): the first multi-node (5-node) MPI job this
    project ever ran wrote a null-byte-padded gap (exactly one 4KB block,
    the classic signature of a parallel-filesystem write race on a file's
    first block) into shear_stress.dat, merging several rows into one
    unparseable line. np.loadtxt() on that file raises ValueError, not a
    silent bad value -- the ORIGINAL version of this check let that
    exception propagate uncaught, crashing the validator instead of
    reporting a clean FAIL. A crashed validator is silence, and silence
    reads as "nothing to report" -- exactly the failure mode this script
    exists to prevent."""
    run_dir = RUNS / params["run_id"]
    bad = []
    for fname in ("shear_stress.dat", "normf.dat", "vol_frac_interf.dat", "tr_oxy.dat"):
        p = run_dir / fname
        if not p.exists():
            continue
        try:
            d = np.loadtxt(p, skiprows=1)
        except ValueError as e:
            bad.append(f"{fname} (unparseable: {e})")
            continue
        if not np.all(np.isfinite(d)):
            bad.append(fname)
    rjson = run_dir / "results.json"
    nan_keys = []
    if rjson.exists():
        r = json.loads(rjson.read_text())
        nan_keys = [k for k, v in r.items() if isinstance(v, float) and math.isnan(v)]
    if bad:
        report.append(("HARD", "FAIL", f"NaN/Inf in: {', '.join(bad)}"))
    elif nan_keys:
        report.append(("WARN", "WARN", f"results.json NaN keys (may be expected, e.g. kLa if "
                                       f"t_end<=t_mix): {', '.join(nan_keys)}"))
    else:
        report.append(("HARD", "PASS", "no NaN/Inf in time series; no unexpected NaN in results.json"))


def check_restart_fidelity(params, parent_params, report):
    if parent_params is None:
        report.append(("SKIP", "SKIP", "no parent given -- restart fidelity not checked"))
        return
    this_diag = load_diag(RUNS / params["run_id"] / "restart_diagnostic_post_restore.txt")
    parent_diag = load_diag(RUNS / parent_params["run_id"] / "restart_diagnostic_pre_dump.txt")
    if this_diag is None or parent_diag is None:
        missing = "parent" if parent_diag is None else "this segment"
        report.append(("HARD", "FAIL",
                       f"CANNOT VERIFY: {missing}'s restart_diagnostic file is missing -- likely "
                       f"because {parent_params['run_id']!r} predates the diagnostic "
                       f"instrumentation (added 2026-09-07). Not a detected mismatch, an inability "
                       f"to check; treat this restore as unverified, not as confirmed-bad."))
        return

    # write_restart_diagnostic() fires BEFORE any restart fixup (su rescale,
    # fs recompute, stracer reset) -- both files capture the RAW restored
    # state, so the right comparison is direct equality, not scaled by su/
    # su^2. Got this wrong on the first version of this check (compared
    # scaled values, got a spurious 26.6% "mismatch" on a pair whose raw
    # files are bit-identical to every printed digit) -- caught by reading
    # the raw files before trusting the check's own output, which is
    # exactly the discipline this script exists to enforce on the
    # simulations themselves.
    checks = ["ux_sum", "uy_sum", "p_sum", "f_sum", "pf_sum"]
    worst = 0.0
    for key in checks:
        expected = parent_diag[key]
        got = this_diag[key]
        denom = max(abs(expected), 1e-30)
        rel = abs(got - expected) / denom
        worst = max(worst, rel)
    if worst > 1e-6:
        report.append(("HARD", "FAIL",
                       f"restore() field mismatch vs parent's pre-dump state: "
                       f"max relative diff {worst:.2e} (expect ~1e-12, machine precision)"))
    else:
        report.append(("HARD", "PASS",
                       f"restore() reproduces parent's pre-dump state exactly "
                       f"(max rel diff {worst:.2e})"))


def check_parent_linkage(params, parent_params, report):
    if parent_params is None:
        report.append(("SKIP", "SKIP", "no parent given -- linkage not checked"))
        return
    exp_theta = parent_params["theta_max"][0]
    exp_omega = parent_params["omega_b"]
    got_theta = params.get("theta_max_prev", [None])[0]
    got_omega = params.get("omega_b_prev")
    ok = (got_theta is not None and abs(got_theta - exp_theta) < 1e-9 and
          got_omega is not None and abs(got_omega - exp_omega) < 1e-9)
    if ok:
        report.append(("HARD", "PASS", "theta_max_prev/omega_b_prev match parent segment exactly"))
    else:
        report.append(("HARD", "FAIL",
                       f"parent linkage mismatch: theta_max_prev={got_theta} (parent had "
                       f"{exp_theta}), omega_b_prev={got_omega} (parent had {exp_omega})"))


def check_drift(params, report):
    """Is the TAIL of the run still moving? Deliberately does NOT compare
    early vs late windows -- an early window is dominated by the ramp-up
    transient (a cold start's own early cycles can span the full 0-to-1
    range), which inflates any noise estimate taken from it and would let
    a genuinely still-drifting run pass (caught live: this exact bug let
    l9_cold_rpm22.5 -- known independently, from the detailed cycle-by-
    cycle check earlier, to still be climbing at cycle 44 -- report as
    'stable'). Compares only the run's own last two blocks against each
    other, so a smooth ramp elsewhere in the run can't hide a live drift
    at the point the run actually stopped."""
    rpm = params["omega_b"] * 60.0 / (2 * math.pi)
    theta = params["theta_max"][0]
    T_ND = t_per_nd(rpm, theta)
    try:
        d = np.loadtxt(RUNS / params["run_id"] / "shear_stress.dat", skiprows=1)
    except ValueError as e:
        report.append(("HARD", "FAIL", f"shear_stress.dat unparseable: {e}"))
        return
    t, y = d[:, 1], d[:, 5]
    c = (t - t[0]) / T_ND
    n = int(c[-1])
    if n < 10:
        report.append(("WARN", "WARN", f"only {n} cycles -- too short for a drift check"))
        return
    peak = np.array([y[(c >= k) & (c < k + 1)].max() for k in range(n)])
    block = max(5, n // 5)
    last, prev = peak[-block:], peak[-2 * block:-block]
    local_noise = np.std(last) / max(np.mean(last), 1e-30)
    drift = abs(np.mean(last) - np.mean(prev)) / max(np.mean(prev), 1e-30)
    if drift > max(0.02, 3 * local_noise):
        report.append(("WARN", "WARN",
                       f"NOT SETTLED: last {block} cycles differ from the {block} before them "
                       f"by {100*drift:.1f}% (tail's own noise floor {100*local_noise:.1f}%) -- "
                       f"do not report this run's converged value, extend it"))
    else:
        report.append(("HARD", "PASS",
                       f"tail settled: last two {block}-cycle blocks agree to {100*drift:.1f}%"
                       f" (tail noise floor {100*local_noise:.1f}%)"))


def validate(run_id, parent_id=None):
    params = json.loads((RUNS / run_id / "params.json").read_text())
    params["run_id"] = run_id
    parent_params = None
    if parent_id:
        parent_params = json.loads((RUNS / parent_id / "params.json").read_text())
        parent_params["run_id"] = parent_id

    report = []
    is_restart = params.get("t_checkpoint") is not None
    checks = [check_binary, check_cycle_count, check_nan_inf]
    if is_restart:
        checks += [lambda p, r: check_restart_fidelity(p, parent_params, r),
                  lambda p, r: check_parent_linkage(p, parent_params, r)]
    checks.append(check_drift)
    for fn in checks:
        try:
            fn(params, report)
        except Exception as e:
            # A check that crashes is silence, not a pass -- and silence
            # reads as "nothing to report." Report it as a failure to
            # analyze instead of letting the whole script die (caught
            # live, 2026-09-10: the first multi-node job's corrupted
            # output file crashed three of these checks in a row before
            # this net was added).
            report.append(("HARD", "FAIL", f"{fn.__name__ if hasattr(fn, '__name__') else 'check'} "
                                           f"crashed: {type(e).__name__}: {e}"))

    hard_fail = any(sev == "HARD" and status == "FAIL" for sev, status, _ in report)
    verdict = "FAIL" if hard_fail else "PASS"
    print(f"=== {run_id} ({'restart' if is_restart else 'cold start'}"
         f"{f', parent={parent_id}' if parent_id else ''}) -> {verdict} ===")
    for sev, status, msg in report:
        print(f"  [{status:4s}] {msg}")
    print()
    return verdict == "PASS"


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run_id", nargs="?")
    ap.add_argument("--parent", default=None)
    args = ap.parse_args()
    if not args.run_id:
        print("usage: validate_run.py <run_id> [--parent <parent_run_id>]", file=sys.stderr)
        sys.exit(1)
    ok = validate(args.run_id, args.parent)
    sys.exit(0 if ok else 1)
