"""Regression guards: MPI parallelism and checkpoint-restart must not change
physics, for OUR OWN fork, relative to OUR OWN serial/uninterrupted baseline.

Why these tests exist
----------------------
2026-07-30 investigation (diary.md) found that our fork's serial baseline
disagrees with a near-literal reproduction of Kim et al.'s own upstream code
by ~1.74x (5.42 vs 9.44 in ux_rms/U_b, same condition/resolution) -- and
manually confirmed MPI and checkpoint-restart are NOT the cause (both agree
with the fork's own serial baseline to within ~2%). These tests turn that
one-off manual investigation into a standing regression guard, run
automatically by CI: if a future change to the MPI or checkpoint-restart
code path silently breaks physics, `medium-tests` catches it, without
being confused by (and without ever asserting anything about) how our
fork compares to Kim's own code.

CI status: marked `medium` and DOES run in GitHub Actions' `medium-tests`
job (ci.yml installs openmpi via apt and runs `make build-mpi` before this
suite -- there is no SLURM/module-system dependency; `make build-mpi` was
fixed to skip `module load` when `mpicc` is already on PATH). The
MPI-dependent tests here (`test_mpi_matches_serial`,
`test_combined_mpi_checkpoint_vs_serial`) skip gracefully via the
`mpi_binary` fixture if `make build-mpi` fails for any reason;
`test_checkpoint_matches_uninterrupted` needs neither MPI nor SLURM at
all -- it's a plain serial binary run twice.

Baseline policy (explicit user decision, 2026-07-30): the MANDATORY,
pass/fail assertions in this file compare our fork against ITSELF (serial
vs MPI, uninterrupted vs checkpoint-restarted) -- never against Kim's
upstream code. A separate, non-fatal WARNING-only comparison against a
vendored copy of Kim's own driver code lives in
test_kim_upstream_comparison.py; see that file's docstring for why it never
raises.

All tests here are real CFD runs at fidelity 5 (32x32 cells) or fidelity 6
where a checkpoint fixture is reused from another test.
"""
from __future__ import annotations

import json
import math
import shutil
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parents[2]))
from tests.conftest import CANONICAL_PARAMS, load_normf, run_bioreactor
from tests.verification.test_grid_convergence import _t_ramp_nd, _mean_post_ramp_vel_rms

PROJECT_ROOT = Path(__file__).parents[2]

# Cheap-but-real condition for MPI/checkpoint parity checks: same physical
# condition as CANONICAL_PARAMS, fidelity dropped to 5 (32x32) to keep these
# fast enough to run routinely on a compute node.
_BASE_PARAMS = {**CANONICAL_PARAMS, "fidelity": 5, "t_end": 12.0}
_TIMEOUT = 1800  # seconds; fidelity=5 is much cheaper than the L5 grid-convergence test's 7200s

_VEL_RTOL_MPI  = 0.05   # MPI vs serial: manually measured ~0.2% -- 5% leaves large margin
_VEL_RTOL_CKPT = 0.10   # checkpoint vs uninterrupted: manually measured ~1.6% -- 10% margin


@pytest.fixture(scope="module")
def mpi_binary():
    """Build BioReactor-mpi via the project Makefile. Skips if openmpi unavailable."""
    build = subprocess.run(
        ["make", "build-mpi"], cwd=PROJECT_ROOT, capture_output=True, text=True,
    )
    binary = PROJECT_ROOT / "build" / "BioReactor-mpi"
    if build.returncode != 0 or not binary.exists():
        pytest.skip(f"build-mpi failed (openmpi module likely unavailable):\n{build.stderr}")
    return binary


def _run_mpi(binary: Path, params: dict, tmp_path: Path, nranks: int, timeout: int) -> Path:
    run_dir = tmp_path / params["run_id"]
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "params.json").write_text(json.dumps(params))
    try:
        subprocess.run(
            ["mpirun", "-np", str(nranks), "--oversubscribe", str(binary.resolve()), "params.json"],
            cwd=run_dir, capture_output=True, text=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        pass
    return run_dir


def _vel_rms(run_dir: Path, params: dict) -> float:
    normf_path = run_dir / "normf.dat"
    if not normf_path.exists():
        pytest.fail(f"normf.dat not written in {run_dir} -- simulation likely crashed")
    return _mean_post_ramp_vel_rms(load_normf(run_dir), _t_ramp_nd(params))


@pytest.mark.medium
def test_mpi_matches_serial(mpi_binary, tmp_path):
    """(a) MPI must reproduce our own serial run's post-ramp velocity RMS.

    Mandatory assertion, baseline = our own fork's serial run (not Kim's
    upstream code -- see module docstring). A failure here means MPI domain
    decomposition (not Kim-vs-fork physics differences) broke something.
    """
    serial_params = {**_BASE_PARAMS, "run_id": "mpi_parity_serial"}
    mpi_params    = {**_BASE_PARAMS, "run_id": "mpi_parity_mpi"}

    serial_dir = run_bioreactor(serial_params, tmp_path, timeout=_TIMEOUT)
    mpi_dir    = _run_mpi(mpi_binary, mpi_params, tmp_path, nranks=4, timeout=_TIMEOUT)

    vel_serial = _vel_rms(serial_dir, serial_params)
    vel_mpi    = _vel_rms(mpi_dir, mpi_params)

    rel_err = abs(vel_mpi - vel_serial) / (vel_serial + 1e-30)
    assert rel_err < _VEL_RTOL_MPI, (
        f"MPI velocity RMS diverges from serial: serial={vel_serial:.5f}, mpi={vel_mpi:.5f}, "
        f"relative error={rel_err:.2%} (threshold {_VEL_RTOL_MPI:.0%}). "
        "This means MPI domain decomposition is changing physics, not just floating-point "
        "reduction order -- investigate before trusting any MPI production run."
    )


@pytest.mark.medium
def test_checkpoint_matches_uninterrupted(tmp_path):
    """(b) A checkpoint-restart (same condition across the boundary) must
    reproduce an uninterrupted run's post-ramp velocity RMS.

    Mandatory assertion, baseline = our own fork's uninterrupted serial run.
    `*_prev` params are set EQUAL to the current condition so the fork's own
    smooth-step continuity ramp ((1-alpha)*prev + alpha*current) is a no-op
    -- this isolates pure checkpoint-restart mechanics from any condition
    change across the boundary.
    """
    uninterrupted_params = {**_BASE_PARAMS, "run_id": "ckpt_parity_uninterrupted"}
    uninterrupted_dir = run_bioreactor(uninterrupted_params, tmp_path, timeout=_TIMEOUT)
    vel_uninterrupted = _vel_rms(uninterrupted_dir, uninterrupted_params)

    seg1_params = {**_BASE_PARAMS, "run_id": "ckpt_parity_seg1", "t_end": 6.0}
    seg1_dir = run_bioreactor(seg1_params, tmp_path, timeout=_TIMEOUT)
    dump_path = seg1_dir / "checkpoint.dump"
    if not dump_path.exists():
        pytest.fail(f"checkpoint.dump not written in {seg1_dir} -- seg1 run likely crashed")

    # actual checkpoint time (rounded to a period boundary by the C code) --
    # read it back from seg1's own normf.dat rather than assuming t_end exactly.
    seg1_t_final = float(load_normf(seg1_dir)[-1, 1])

    seg2_params = {
        **_BASE_PARAMS, "run_id": "ckpt_parity_seg2", "t_end": 6.0,
        "t_checkpoint": seg1_t_final,
        "omega_b_prev": _BASE_PARAMS["omega_b"],
        "theta_max_prev": _BASE_PARAMS["theta_max"],
    }
    seg2_dir = tmp_path / seg2_params["run_id"]
    seg2_dir.mkdir(parents=True, exist_ok=True)
    (seg2_dir / "params.json").write_text(json.dumps(seg2_params))
    shutil.copy(dump_path, seg2_dir / "checkpoint.dump")
    binary = PROJECT_ROOT / "build" / "BioReactor"
    try:
        subprocess.run(
            [str(binary.resolve()), "params.json", "checkpoint.dump"],
            cwd=seg2_dir, capture_output=True, text=True, timeout=_TIMEOUT,
        )
    except subprocess.TimeoutExpired:
        pass

    vel_ckpt = _vel_rms(seg2_dir, seg2_params)

    rel_err = abs(vel_ckpt - vel_uninterrupted) / (vel_uninterrupted + 1e-30)
    assert rel_err < _VEL_RTOL_CKPT, (
        f"Checkpoint-restarted velocity RMS diverges from uninterrupted run: "
        f"uninterrupted={vel_uninterrupted:.5f}, checkpoint-restarted={vel_ckpt:.5f}, "
        f"relative error={rel_err:.2%} (threshold {_VEL_RTOL_CKPT:.0%}). "
        "This means the checkpoint dump/restore path itself is altering the physical state, "
        "not just contributing benign floating-point noise -- see the project's own history "
        "of checkpoint-restart correctness bugs (commit 19c3a31) before trusting a chained run."
    )


@pytest.mark.medium
def test_checkpoint_theta_change_matches_independent_reference(tmp_path):
    """Regression test for the 2026-09-07 cold-start-beats-warm-start
    anomaly (diary.md): a checkpoint-restart that changes theta_max mid-run
    should converge to the SAME quasi-steady velocity field as an
    INDEPENDENT run started fresh at the target theta_max -- checkpointing
    exists to skip a redundant transient, not to reach a different state.

    Found instead (L6, 32.5rpm, theta 7->4, `experiments/
    dtheta_monotonicity/settling_vs_dtheta_th4target.png`): the checkpoint-
    restarted run does NOT slowly settle toward the independent reference,
    it PLATEAUS at a persistent 14-23% offset from it, for every nonzero
    theta step tested (even Delta_theta_max=0.1deg) -- while the pure
    identity restart (Delta_theta_max=0, see test_checkpoint_matches_
    uninterrupted above) matches its own reference almost exactly. Two
    hypotheses for the mechanism have already been tested and FALSIFIED
    (see diary.md): (1) a directional "excess energy" story -- ruled out,
    reverse-direction warm starts are equally pathological; (2) a missing
    dAk/dt chain-rule term in the ramp's Th_d/Th_2d formulas -- fixed
    (correct regardless), but the offset was unchanged (16.6% -> 15.3%) on
    the smallest anomalous case. Root cause is still open.

    THIS TEST IS EXPECTED TO FAIL until that anomaly is root-caused and
    fixed -- it exists to formally document the bug at CI scale and serve
    as the eventual fix's acceptance criterion. A green run of this test
    (without a corresponding diary.md entry explaining the fix) should be
    treated as suspicious, not celebrated.
    """
    theta_from, theta_to = 7.0, 4.0

    source_params = {**_BASE_PARAMS, "run_id": "theta_ckpt_source",
                      "theta_max": [theta_from, 0.0, 0.0]}
    source_dir = run_bioreactor(source_params, tmp_path, timeout=_TIMEOUT)
    source_t_final = float(load_normf(source_dir)[-1, 1])
    dump_path = source_dir / "checkpoint.dump"
    if not dump_path.exists():
        pytest.fail(f"checkpoint.dump not written in {source_dir} -- source run likely crashed")

    reference_params = {**_BASE_PARAMS, "run_id": "theta_ckpt_reference",
                         "theta_max": [theta_to, 0.0, 0.0]}
    reference_dir = run_bioreactor(reference_params, tmp_path, timeout=_TIMEOUT)
    vel_reference = _vel_rms(reference_dir, reference_params)

    transition_params = {
        **_BASE_PARAMS, "run_id": "theta_ckpt_transition",
        "theta_max": [theta_to, 0.0, 0.0], "t_end": 12.0,
        "t_checkpoint": source_t_final,
        "omega_b_prev": _BASE_PARAMS["omega_b"],
        "theta_max_prev": [theta_from, 0.0, 0.0],
    }
    transition_dir = tmp_path / transition_params["run_id"]
    transition_dir.mkdir(parents=True, exist_ok=True)
    (transition_dir / "params.json").write_text(json.dumps(transition_params))
    shutil.copy(dump_path, transition_dir / "checkpoint.dump")
    binary = PROJECT_ROOT / "build" / "BioReactor"
    try:
        subprocess.run(
            [str(binary.resolve()), "params.json", "checkpoint.dump"],
            cwd=transition_dir, capture_output=True, text=True, timeout=_TIMEOUT,
        )
    except subprocess.TimeoutExpired:
        pass

    # Ramp for the transition segment starts AT t_checkpoint (absolute sim
    # time never resets on restart), not at t=0 -- _t_ramp_nd() alone would
    # under-window here (it returns a 3-cycle DURATION, meant to be added to
    # the segment's own start time, not compared to it directly).
    normf_transition = load_normf(transition_dir)
    t_ramp_abs = source_t_final + _t_ramp_nd(transition_params)
    vel_transition = _mean_post_ramp_vel_rms(normf_transition, t_ramp_abs)

    assert not math.isnan(vel_reference) and not math.isnan(vel_transition), (
        f"vel_rms is NaN (reference={vel_reference}, transition={vel_transition}) -- "
        "one of the runs produced too little post-ramp data."
    )
    rel_err = abs(vel_transition - vel_reference) / (vel_reference + 1e-30)
    assert rel_err < _VEL_RTOL_CKPT, (
        f"Checkpoint-restarted velocity RMS after a theta_max change diverges from an "
        f"INDEPENDENT reference run at the same target condition: "
        f"reference={vel_reference:.5f}, transition={vel_transition:.5f}, "
        f"relative error={rel_err:.2%} (threshold {_VEL_RTOL_CKPT:.0%}). "
        "See diary.md 2026-09-07: this is the cold-start-beats-warm-start anomaly, "
        "root cause still open as of this test's creation."
    )


@pytest.mark.medium
def test_combined_mpi_checkpoint_vs_serial(mpi_binary, tmp_path):
    """(c) MPI + checkpoint-restart TOGETHER (the actual production
    configuration) must reproduce a plain serial/uninterrupted run's
    velocity field and stress field. kLa is reported but not asserted on --
    kLa is known to be highly sensitive to small perturbations at low
    fidelity even between nominally-identical runs (see
    test_grid_convergence.py's own kLa caveat); asserting a tight tolerance
    on it here would make this test flaky for reasons unrelated to whether
    MPI+checkpoint corrupts anything.
    """
    from scripts import postprocess

    combo_params = {
        **_BASE_PARAMS, "run_id": "combo_serial_baseline",
        "n_mix_cycles": 5, "t_end": 10.0,
    }
    baseline_dir = run_bioreactor(combo_params, tmp_path, timeout=_TIMEOUT)
    baseline_kpis = postprocess.main(str(baseline_dir), params=combo_params)

    # MPI + checkpoint together: split the same condition into two MPI segments.
    seg1_params = {**combo_params, "run_id": "combo_mpi_ckpt_seg1", "t_end": 5.0}
    seg1_dir = _run_mpi(mpi_binary, seg1_params, tmp_path, nranks=4, timeout=_TIMEOUT)
    dump_path = seg1_dir / "checkpoint.dump"
    if not dump_path.exists():
        pytest.fail(f"checkpoint.dump not written in {seg1_dir}")
    seg1_t_final = float(load_normf(seg1_dir)[-1, 1])

    seg2_params = {
        **combo_params, "run_id": "combo_mpi_ckpt_seg2", "t_end": 5.0,
        "t_checkpoint": seg1_t_final,
        "omega_b_prev": combo_params["omega_b"],
        "theta_max_prev": combo_params["theta_max"],
    }
    seg2_dir = tmp_path / seg2_params["run_id"]
    seg2_dir.mkdir(parents=True, exist_ok=True)
    (seg2_dir / "params.json").write_text(json.dumps(seg2_params))
    shutil.copy(dump_path, seg2_dir / "checkpoint.dump")
    try:
        subprocess.run(
            ["mpirun", "-np", "4", "--oversubscribe", str(mpi_binary.resolve()),
             "params.json", "checkpoint.dump"],
            cwd=seg2_dir, capture_output=True, text=True, timeout=_TIMEOUT,
        )
    except subprocess.TimeoutExpired:
        pass
    combo_kpis = postprocess.main(str(seg2_dir), params=seg2_params)

    # ── velocity field (mandatory) ──────────────────────────────────────────
    vel_base, vel_combo = baseline_kpis["vel_rms_qss"], combo_kpis["vel_rms_qss"]
    assert not math.isnan(vel_base) and not math.isnan(vel_combo), (
        f"vel_rms_qss is NaN (baseline={vel_base}, combo={vel_combo}) -- "
        "one of the runs produced too little post-ramp/pre-injection data."
    )
    vel_rel_err = abs(vel_combo - vel_base) / (vel_base + 1e-30)
    assert vel_rel_err < 0.15, (
        f"Combined MPI+checkpoint velocity field diverges from serial baseline: "
        f"baseline={vel_base:.5f}, combo={vel_combo:.5f}, rel_err={vel_rel_err:.2%} "
        "(threshold 15%)."
    )

    # ── stress field (mandatory) ────────────────────────────────────────────
    # Use tau_mean_max (max-over-time of the SPATIALLY-AVERAGED stress), not
    # tau_100_max (absolute max over all space AND time): tau_100_max is a
    # pure extreme-value statistic from a short, coarsely-sampled run and is
    # inherently much noisier -- whether the discrete output cadence happens
    # to land on the true instantaneous peak differs between a continuous
    # 10-unit run and a 5+5 checkpoint-split run even with identical physics.
    # tau_mean_max is far smoother and a fairer test of "did MPI/checkpoint
    # change the physics" rather than "did sampling luck change".
    tau_base, tau_combo = baseline_kpis["tau_mean_max"], combo_kpis["tau_mean_max"]
    if not (math.isnan(tau_base) or math.isnan(tau_combo)):
        tau_rel_err = abs(tau_combo - tau_base) / (tau_base + 1e-30)
        assert tau_rel_err < 0.20, (
            f"Combined MPI+checkpoint stress field (tau_mean_max) diverges from serial "
            f"baseline: baseline={tau_base:.5f}, combo={tau_combo:.5f}, "
            f"rel_err={tau_rel_err:.2%} (threshold 20%)."
        )
    else:
        warnings.warn(
            f"tau_mean_max is NaN for baseline={tau_base} or combo={tau_combo} -- "
            "skipping the stress-field assertion (insufficient shear_stress.dat data)."
        )

    # tau_100_max (extreme-value statistic) reported for visibility, never asserted --
    # see comment above on why it's expected to be noisier than tau_mean_max.
    tau100_base, tau100_combo = baseline_kpis.get("tau_100_max"), combo_kpis.get("tau_100_max")
    if tau100_base and tau100_combo and not (math.isnan(tau100_base) or math.isnan(tau100_combo)):
        tau100_rel_err = abs(tau100_combo - tau100_base) / (tau100_base + 1e-30)
        if tau100_rel_err > 0.30:
            warnings.warn(
                f"tau_100_max differs by {tau100_rel_err:.1%} between combined MPI+checkpoint "
                f"({tau100_combo:.5f}) and serial baseline ({tau100_base:.5f}). Not asserted on "
                "-- this is an extreme-value statistic, expected to be sensitive to sampling "
                "cadence even with identical physics (see tau_mean_max above for the real check)."
            )

    # ── kLa (reported, non-fatal) ───────────────────────────────────────────
    kla_base, kla_combo = baseline_kpis.get("kLa_25"), combo_kpis.get("kLa_25")
    if kla_base and kla_combo and not (math.isnan(kla_base) or math.isnan(kla_combo)):
        kla_rel_err = abs(kla_combo - kla_base) / (kla_base + 1e-30)
        if kla_rel_err > 0.5:
            warnings.warn(
                f"kLa_25 differs by {kla_rel_err:.1%} between combined MPI+checkpoint "
                f"({kla_combo:.4f}) and serial baseline ({kla_base:.4f}) at fidelity=5, "
                "n_mix_cycles=5. Not asserted on -- kLa is known to be highly sensitive to "
                "small perturbations at low fidelity (see test_grid_convergence.py) -- but "
                "worth a human glance if this warning appears on a change to the MPI/"
                "checkpoint code path specifically."
            )
    else:
        warnings.warn(
            f"kLa_25 unavailable (baseline={kla_base}, combo={kla_combo}) -- oxygen "
            "likely never reached 25% saturation in this short test window."
        )
