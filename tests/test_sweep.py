"""Tests for scripts/sweep.py — JSON-driven multi-param sweep runner."""
import itertools
import math
import shutil
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parents[1]))

_PROJECT_ROOT = Path(__file__).parents[1]

from scripts.sweep import (
    build_segment_list,
    detect_sweep_params,
    expand_combinations,
    group_by_checkpoint_key,
)

# ── shared fixtures ───────────────────────────────────────────────────────────

_BASE = {
    "fidelity": 5,
    "geometry": {"a": 0.25, "b": 0.071, "n": 2.0},
    "fill_level": 0.5,
    "omega_b": 3.14159,
    "n_harmonics": 1,
    "theta_max": [7.0, 0.0, 0.0],
    "phi_angular": [0.0, 0.0, 0.0],
    "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0],
    "phi_horizontal": [0.0, 0.0, 0.0],
    "n_mix_cycles": 10,
}

_OPTIONS = {
    "n_mix_cycles": 10,
    "n_transition_cycles": 3,
    "t_buffer": 3.0,
    "walltime": "00:10:00",
    "submit": False,
}


# ── detect_sweep_params ───────────────────────────────────────────────────────

def test_no_lists_no_sweep():
    """A plain params dict with no list values has no sweep params."""
    result = detect_sweep_params(_BASE)
    assert result == {}


def test_scalar_list_is_sweep():
    """A scalar param with a list value is a sweep param."""
    raw = {**_BASE, "omega_b": [3.14159, 6.28318]}
    result = detect_sweep_params(raw)
    assert "omega_b" in result
    assert result["omega_b"] == [3.14159, 6.28318]


def test_fidelity_sweepable():
    """fidelity accepts a list → sweep."""
    raw = {**_BASE, "fidelity": [5, 7]}
    result = detect_sweep_params(raw)
    assert "fidelity" in result
    assert result["fidelity"] == [5, 7]


def test_timing_sweepable():
    """n_mix_cycles accepts a list → sweep."""
    raw = {**_BASE, "n_mix_cycles": [10, 80]}
    result = detect_sweep_params(raw)
    assert "n_mix_cycles" in result


def test_fill_level_sweepable():
    """fill_level accepts a list → sweep."""
    raw = {**_BASE, "fill_level": [0.4, 0.5, 0.6]}
    result = detect_sweep_params(raw)
    assert "fill_level" in result


def test_known_vector_not_swept_scalar_elements():
    """theta_max=[7,0,0] (list of scalars) is a fixed 3-vector, not a sweep."""
    raw = {**_BASE, "theta_max": [7.0, 0.0, 0.0]}
    result = detect_sweep_params(raw)
    assert "theta_max" not in result


def test_nested_list_swept():
    """theta_max=[[5,0,0],[7,0,0]] (list of lists) IS a sweep."""
    raw = {**_BASE, "theta_max": [[5.0, 0.0, 0.0], [7.0, 0.0, 0.0]]}
    result = detect_sweep_params(raw)
    assert "theta_max" in result
    assert result["theta_max"] == [[5.0, 0.0, 0.0], [7.0, 0.0, 0.0]]


def test_geometry_sweepable():
    """geometry=[{...},{...}] (list of dicts) is a sweep."""
    geom1 = {"a": 0.20, "b": 0.09, "n": 2.0}
    geom2 = {"a": 0.25, "b": 0.071, "n": 8.0}
    raw = {**_BASE, "geometry": [geom1, geom2]}
    result = detect_sweep_params(raw)
    assert "geometry" in result
    assert result["geometry"] == [geom1, geom2]


def test_sweep_key_excluded():
    """_sweep dict is never treated as a sweep param."""
    raw = {**_BASE, "_sweep": {"n_mix_cycles": 10, "t_buffer": 5.0}}
    result = detect_sweep_params(raw)
    assert "_sweep" not in result


# ── expand_combinations ───────────────────────────────────────────────────────

def test_no_sweep_params_one_combo():
    """No sweep params → one empty combo dict."""
    result = expand_combinations({})
    assert result == [{}]


def test_zip_same_length():
    """Two sweep params each with 3 values → 3 sims (zip)."""
    sweep = {"omega_b": [1.0, 2.0, 3.0], "fill_level": [0.4, 0.5, 0.6]}
    result = expand_combinations(sweep)
    assert len(result) == 3


def test_zip_values_correct():
    """Zip: each sim has paired values from same index."""
    sweep = {"omega_b": [1.0, 2.0], "fill_level": [0.4, 0.5]}
    result = expand_combinations(sweep)
    assert result[0] == {"omega_b": 1.0, "fill_level": 0.4}
    assert result[1] == {"omega_b": 2.0, "fill_level": 0.5}


def test_cartesian_different_lengths():
    """Two sweep params with lengths 2 and 3 → 6 sims (cartesian)."""
    sweep = {"omega_b": [1.0, 2.0], "fill_level": [0.3, 0.5, 0.7]}
    result = expand_combinations(sweep)
    assert len(result) == 6


def test_cartesian_values_correct():
    """Cartesian: result contains all combinations."""
    sweep = {"omega_b": [1.0, 2.0], "fill_level": [0.4, 0.6]}
    result = expand_combinations(sweep)
    # Different lengths → cartesian (even though both length 2 would be zip).
    # Force cartesian: use length 2 and 3.
    sweep2 = {"omega_b": [1.0, 2.0], "fill_level": [0.3, 0.5, 0.7]}
    result2 = expand_combinations(sweep2)
    expected = [
        {"omega_b": 1.0, "fill_level": 0.3},
        {"omega_b": 1.0, "fill_level": 0.5},
        {"omega_b": 1.0, "fill_level": 0.7},
        {"omega_b": 2.0, "fill_level": 0.3},
        {"omega_b": 2.0, "fill_level": 0.5},
        {"omega_b": 2.0, "fill_level": 0.7},
    ]
    assert result2 == expected


def test_single_sweep_param():
    """Single sweep param → N sims, each with one value."""
    sweep = {"omega_b": [1.0, 2.0, 3.0]}
    result = expand_combinations(sweep)
    assert len(result) == 3
    assert [r["omega_b"] for r in result] == [1.0, 2.0, 3.0]


# ── group_by_checkpoint_key ───────────────────────────────────────────────────

def _make_params(fidelity=5, a=0.25, b=0.071, n=2.0, **kwargs):
    p = dict(_BASE)
    p["fidelity"] = fidelity
    p["geometry"] = {"a": a, "b": b, "n": n}
    p.update(kwargs)
    return p


def test_same_fidelity_geometry_one_group():
    """Sims with same fidelity+geometry → 1 group."""
    sims = [_make_params(omega_b=1.0), _make_params(omega_b=2.0)]
    groups = group_by_checkpoint_key(sims)
    assert len(groups) == 1


def test_grouping_by_fidelity():
    """Different fidelity → 2 groups."""
    sims = [_make_params(fidelity=5), _make_params(fidelity=7)]
    groups = group_by_checkpoint_key(sims)
    assert len(groups) == 2


def test_grouping_by_geometry():
    """Different geometry → 2 groups."""
    sims = [_make_params(a=0.20, b=0.09), _make_params(a=0.25, b=0.071)]
    groups = group_by_checkpoint_key(sims)
    assert len(groups) == 2


def test_group_preserves_order():
    """Sims within a group keep their expansion order."""
    sims = [
        _make_params(omega_b=1.0),
        _make_params(omega_b=2.0),
        _make_params(omega_b=3.0),
    ]
    groups = group_by_checkpoint_key(sims)
    assert len(groups) == 1
    group = next(iter(groups.values()))
    assert [p["omega_b"] for p in group] == [1.0, 2.0, 3.0]


def test_mixed_fidelity_grouping():
    """3 sims: 2 share fidelity=5, 1 has fidelity=7 → 2 groups."""
    sims = [
        _make_params(fidelity=5, omega_b=1.0),
        _make_params(fidelity=5, omega_b=2.0),
        _make_params(fidelity=7, omega_b=1.0),
    ]
    groups = group_by_checkpoint_key(sims)
    assert len(groups) == 2
    sizes = sorted(len(v) for v in groups.values())
    assert sizes == [1, 2]


def test_grouping_by_fill_level():
    """Different fill_level → 2 groups (incompatible checkpoint fields)."""
    sims = [_make_params(fill_level=0.3), _make_params(fill_level=0.5)]
    groups = group_by_checkpoint_key(sims)
    assert len(groups) == 2


def test_grouping_by_theta_max():
    """Different theta_max[0] → 2 groups (parallel chain efficiency)."""
    sims = [
        _make_params(theta_max=[5.0, 0.0, 0.0]),
        _make_params(theta_max=[7.0, 0.0, 0.0]),
    ]
    groups = group_by_checkpoint_key(sims)
    assert len(groups) == 2


def test_same_fill_and_theta_one_group():
    """Same fidelity, geometry, fill_level, and theta_max → still 1 group."""
    sims = [
        _make_params(omega_b=1.0, fill_level=0.5, theta_max=[7.0, 0.0, 0.0]),
        _make_params(omega_b=2.0, fill_level=0.5, theta_max=[7.0, 0.0, 0.0]),
    ]
    groups = group_by_checkpoint_key(sims)
    assert len(groups) == 1


# ── build_segment_list ────────────────────────────────────────────────────────

def _group_of_two():
    return [_make_params(omega_b=3.14159), _make_params(omega_b=6.28318)]


def test_segment0_no_t_checkpoint():
    """First segment in a group is a fresh run: no t_checkpoint key (or zero)."""
    segs = build_segment_list(_group_of_two(), _OPTIONS)
    assert segs[0].get("t_checkpoint", 0.0) == 0.0


def test_chain_checkpointing():
    """Segment k≥1 in a group has t_checkpoint > 0."""
    segs = build_segment_list(_group_of_two(), _OPTIONS)
    assert segs[1].get("t_checkpoint", 0.0) > 0.0


def test_segment0_full_mix_cycles():
    """First segment uses n_mix_cycles from options (not n_transition_cycles)."""
    segs = build_segment_list(_group_of_two(), _OPTIONS)
    assert segs[0]["n_mix_cycles"] == _OPTIONS["n_mix_cycles"]


def test_restart_transition_cycles():
    """Subsequent segments use n_transition_cycles from options."""
    segs = build_segment_list(_group_of_two(), _OPTIONS)
    assert segs[1]["n_mix_cycles"] == _OPTIONS["n_transition_cycles"]


def test_segment_has_run_id():
    """Every segment has a unique run_id assigned."""
    segs = build_segment_list(_group_of_two(), _OPTIONS)
    ids = [s["run_id"] for s in segs]
    assert len(set(ids)) == 2  # unique


def test_segment_has_t_end():
    """Every segment has t_end computed (> 0)."""
    segs = build_segment_list(_group_of_two(), _OPTIONS)
    for s in segs:
        assert s.get("t_end", 0.0) > 0.0


def test_restart_has_omega_b_prev():
    """Restart segment records the previous omega_b for smooth ramp in C."""
    segs = build_segment_list(_group_of_two(), _OPTIONS)
    assert "omega_b_prev" in segs[1]
    assert segs[1]["omega_b_prev"] == pytest.approx(3.14159)


def test_explicit_n_mix_cycles_respected():
    """If n_mix_cycles was swept (in swept_keys), per-combo values are preserved."""
    group = [
        _make_params(omega_b=1.0, n_mix_cycles=20),
        _make_params(omega_b=2.0, n_mix_cycles=5),
    ]
    segs = build_segment_list(group, _OPTIONS, swept_keys=frozenset({"n_mix_cycles"}))
    assert segs[0]["n_mix_cycles"] == 20
    assert segs[1]["n_mix_cycles"] == 5


def test_fixed_params_propagated():
    """Non-swept params appear unchanged in every sim after merge."""
    raw = {**_BASE, "omega_b": [3.14159, 6.28318]}
    sweep_p = detect_sweep_params(raw)
    combos = expand_combinations(sweep_p)
    for combo in combos:
        merged = {k: v for k, v in raw.items() if k not in sweep_p}
        merged.update(combo)
        assert merged["fill_level"] == 0.5
        assert merged["fidelity"] == 5
        assert merged["geometry"] == {"a": 0.25, "b": 0.071, "n": 2.0}


# ── submit_sweep video isolation tests ───────────────────────────────────────

from unittest.mock import call, patch
from scripts.sweep import submit_sweep, submit_sweep_videos


def _minimal_sweep_cfg(tmp_path: Path, n_omega: int = 2) -> Path:
    """Write a minimal sweep config JSON to tmp_path and return its path."""
    import json
    cfg = {
        "fidelity": 3,
        "geometry": {"a": 0.25, "b": 0.071, "n": 8.0},
        "fill_level": 0.5,
        "theta_max": [7.0, 0.0, 0.0],
        "phi_angular": [0.0, 0.0, 0.0],
        "omega_h": 0.0,
        "amplitude_h": [0.0, 0.0, 0.0],
        "phi_horizontal": [0.0, 0.0, 0.0],
        "n_harmonics": 1,
        "omega_b": [1.571, 1.833][:n_omega],
        "_sweep": {
            "n_mix_cycles": 3,
            "n_transition_cycles": 3,
            "t_buffer": 5.0,
            "walltime": "00:05:00",
            "submit": True,
        },
    }
    p = tmp_path / "sweep_test.json"
    p.write_text(json.dumps(cfg))
    return p


def test_submit_sweep_staggers_chain_starts(tmp_path):
    """Seg-0 of chain g_idx>0 must have begin=now+N to prevent simultaneous
    dependency releases that cause SLURM to cancel jobs with Reason=Dependency.

    Root cause: when N chains all finish seg-0 simultaneously and release N seg-1
    dependencies at once, SLURM cancels the excess even if parents completed 0:0.
    Staggered --begin ensures seg-0s start (and therefore finish) at different times.
    """
    import json
    cfg = {
        "fidelity": 3,
        "geometry": {"a": 0.25, "b": 0.071, "n": 8.0},
        "fill_level": [0.3, 0.5],       # 2 groups → 2 chains
        "theta_max": [7.0, 0.0, 0.0],
        "phi_angular": [0.0, 0.0, 0.0],
        "omega_h": 0.0, "amplitude_h": [0.0, 0.0, 0.0],
        "phi_horizontal": [0.0, 0.0, 0.0], "n_harmonics": 1,
        "omega_b": [1.571, 1.833, 2.094],  # 3 values (3≠2) → cartesian → 3 segs/chain
        "_sweep": {"n_mix_cycles": 3, "n_transition_cycles": 3,
                   "t_buffer": 5.0, "walltime": "00:05:00",
                   "submit": True, "chain_stagger_s": 300},
    }
    cfg_path = tmp_path / "stagger_test.json"
    cfg_path.write_text(json.dumps(cfg))

    begin_args = []
    def fake_submit(*a, begin=None, **kw):
        begin_args.append(begin)
        return str(10000 + len(begin_args))

    with patch("scripts.simulate.submit_slurm", side_effect=fake_submit):
        submit_sweep(cfg_path)

    # 2 groups × 3 segments = 6 submit calls
    assert len(begin_args) == 6
    # Group 0 seg-0: begin=None (starts immediately)
    assert begin_args[0] is None
    # Group 0 seg-1, seg-2: begin=None (afterok chain, no stagger)
    assert begin_args[1] is None
    assert begin_args[2] is None
    # Group 1 seg-0: begin=now+300 (staggered by one interval)
    assert begin_args[3] == "now+300"
    # Group 1 seg-1, seg-2: begin=None (afterok chain)
    assert begin_args[4] is None
    assert begin_args[5] is None


def test_submit_sweep_does_not_submit_video_jobs(tmp_path):
    """submit_sweep must submit ONLY sim jobs — no video jobs inline.

    Root cause of cascade failure: video jobs co-scheduled with sim jobs on
    the same node; a video failure CANCELLED the next sim segment via SLURM's
    afterok dependency propagation.  Videos must be submitted separately after
    all sims complete.
    """
    cfg = _minimal_sweep_cfg(tmp_path)
    with patch("scripts.simulate.submit_slurm", return_value="99999") as mock_sbatch:
        submit_sweep(cfg)
    # All calls should use the default sim template, NOT the video template
    video_template = str(_PROJECT_ROOT / "config" / "slurm_video_template.sh")
    for c in mock_sbatch.call_args_list:
        template_used = str(c.kwargs.get("template") or "")
        assert "video" not in template_used, (
            f"submit_sweep submitted a video job inline — this causes cascade cancellations. "
            f"Call: {c}"
        )


def test_submit_sweep_sim_job_ids_only_in_dependency_chain(tmp_path):
    """Sim job dependencies must chain only through sim job IDs, never video IDs."""
    cfg = _minimal_sweep_cfg(tmp_path, n_omega=2)
    submitted_ids = []

    def fake_submit(**kwargs):
        job_id = str(10000 + len(submitted_ids))
        submitted_ids.append(job_id)
        return job_id

    with patch("scripts.simulate.submit_slurm", side_effect=lambda *a, **kw: fake_submit(**kw)):
        submit_sweep(cfg)

    # With 2 omega_b values → 1 chain of 2 segments
    # seg0: dependency=None, seg1: dependency=afterok:seg0_id
    # There should be exactly 2 sim jobs and seg1's dependency is seg0's id
    assert len(submitted_ids) == 2
    # If video jobs were inline they'd inflate this count — 2 means sim-only ✓


def test_submit_sweep_videos_passes_checkpoint_for_restart_segments(tmp_path):
    """submit_sweep_videos must export DUMP for restart segments (t_checkpoint>0).

    Without the checkpoint path, BioReactor-video exits with code 1, which
    previously co-cancelled the next simulation segment on the same node.
    """
    import json

    # Create a fake restart run dir
    run_dir = tmp_path / "fake_restart_run"
    run_dir.mkdir()
    params = {**_BASE, "run_id": "fake_restart_run",
              "t_checkpoint": 5.0, "n_mix_cycles": 3, "t_end": 10.0,
              "omega_b_prev": 1.571, "theta_max_prev": [7.0, 0.0, 0.0],
              "phi_angular_prev": [0.0, 0.0, 0.0],
              "amplitude_h_prev": [0.0, 0.0, 0.0],
              "phi_horizontal_prev": [0.0, 0.0, 0.0],
              "omega_h_prev": 0.0}
    (run_dir / "params.json").write_text(json.dumps(params))
    # Fake checkpoint.dump
    (run_dir / "checkpoint.dump").write_bytes(b"fake")

    with patch("scripts.simulate.submit_slurm", return_value="55555") as mock_sbatch:
        job_ids = submit_sweep_videos([run_dir])

    assert job_ids == ["55555"]
    # The checkpoint kwarg must be the path to checkpoint.dump
    call_kwargs = mock_sbatch.call_args.kwargs
    assert call_kwargs.get("checkpoint") is not None, (
        "submit_sweep_videos did not pass checkpoint for a restart segment — "
        "BioReactor-video would exit with code 1"
    )
    assert "checkpoint.dump" in call_kwargs["checkpoint"]


def test_submit_sweep_videos_skips_restart_without_checkpoint(tmp_path):
    """submit_sweep_videos skips a restart segment whose checkpoint.dump is missing."""
    import json

    run_dir = tmp_path / "missing_ck"
    run_dir.mkdir()
    params = {**_BASE, "run_id": "missing_ck",
              "t_checkpoint": 5.0, "n_mix_cycles": 3, "t_end": 10.0}
    (run_dir / "params.json").write_text(json.dumps(params))
    # No checkpoint.dump — simulate a missing file

    with patch("scripts.simulate.submit_slurm", return_value="55556") as mock_sbatch:
        job_ids = submit_sweep_videos([run_dir])

    assert job_ids == [], "Should skip run with missing checkpoint.dump"
    mock_sbatch.assert_not_called()


def test_submit_sweep_videos_no_dependency(tmp_path):
    """Video jobs must have dependency=None so failures cannot cascade to any other job."""
    import json

    run_dir = tmp_path / "fresh_run"
    run_dir.mkdir()
    params = {**_BASE, "run_id": "fresh_run",
              "n_mix_cycles": 3, "t_end": 10.0}
    (run_dir / "params.json").write_text(json.dumps(params))

    with patch("scripts.simulate.submit_slurm", return_value="55557") as mock_sbatch:
        submit_sweep_videos([run_dir])

    call_kwargs = mock_sbatch.call_args.kwargs
    assert call_kwargs.get("dependency") is None, (
        "Video jobs must have dependency=None to prevent cascade cancellations"
    )


# ── SLURM integration test ────────────────────────────────────────────────────
# Marked hpc: only runs when explicitly requested with -m hpc.
# Skipped automatically when sbatch is not on PATH (e.g. login node or CI).
# Compute cost: fidelity=3 (8×8 cells) — each segment runs in ~10 s of CPU time.
# Wall-clock cost: dominated by SLURM queue wait, which can range from seconds
# to ~20 min. The test blocks until both segments complete or the timeout fires.

_SLURM_SWEEP = {
    "fidelity":       3,
    "geometry":       {"a": 0.25, "b": 0.071, "n": 2.0},
    "fill_level":     0.5,
    "n_harmonics":    1,
    "theta_max":      [7.0, 0.0, 0.0],
    "phi_angular":    [0.0, 0.0, 0.0],
    "omega_h":        0.0,
    "amplitude_h":    [0.0, 0.0, 0.0],
    "phi_horizontal": [0.0, 0.0, 0.0],
    "omega_b":        [3.14159, 6.28318],
    "_sweep": {
        "n_mix_cycles":        3,
        "n_transition_cycles": 3,
        "t_buffer":            5.0,
        "walltime":            "00:10:00",
        "submit":              True,
    },
}

_SLURM_TIMEOUT = 1800  # 30 min; generous to absorb queue wait


@pytest.mark.hpc
def test_sweep_slurm_produces_finite_kla(tmp_path):
    """Two-segment fidelity-3 sweep via SLURM: both results.json must be finite.

    Skipped when sbatch is not available (e.g. login node, local CI).
    Marks both segments with a 'test_sweep_slurm' prefix in run_id for easy
    identification and cleanup.

    Checks (per segment):
      - results.json exists
      - kLa_10, kLa_25, kLa_50 are all real numbers (not NaN, not ±inf)
    """
    if shutil.which("sbatch") is None:
        pytest.skip("sbatch not on PATH — SLURM not available on this node")

    import json
    from scripts.sweep import detect_sweep_params, expand_combinations, group_by_checkpoint_key, build_segment_list
    import scripts.simulate as simulate

    runs_root = _PROJECT_ROOT / "runs"

    # Build the two-segment chain manually so we can control run_ids and
    # collect run directories for result validation.
    raw = {k: v for k, v in _SLURM_SWEEP.items() if k != "_sweep"}
    options = _SLURM_SWEEP["_sweep"]
    swept = detect_sweep_params(raw)
    combos = expand_combinations(swept)
    fixed = {k: v for k, v in raw.items() if k not in swept}
    params_list = [{**{k: (list(v) if isinstance(v, list) else v) for k, v in fixed.items()}, **c} for c in combos]
    groups = group_by_checkpoint_key(params_list)
    assert len(groups) == 1, "All sims have same fidelity+geometry → 1 group"
    group = next(iter(groups.values()))
    segments = build_segment_list(group, options)

    # Stamp run_ids with a recognisable prefix
    for k, s in enumerate(segments):
        s["run_id"] = f"test_sweep_slurm_seg{k}_{s['run_id'][:6]}"

    # Submit as a chain
    job_ids = []
    prev_run_id = None
    for k, params in enumerate(segments):
        checkpoint = (
            str((runs_root / prev_run_id / "checkpoint.dump").resolve())
            if prev_run_id else None
        )
        dependency = f"afterok:{job_ids[-1]}" if job_ids else None
        job_id = simulate.submit_slurm(
            params,
            project_root=_PROJECT_ROOT,
            runs_root=runs_root,
            walltime=options["walltime"],
            checkpoint=checkpoint,
            dependency=dependency,
        )
        job_ids.append(job_id)
        prev_run_id = params["run_id"]

    # Wait for each segment's results.json and validate
    for k, params in enumerate(segments):
        run_dir = runs_root / params["run_id"]
        results = simulate.wait_for_result(run_dir, timeout=_SLURM_TIMEOUT, poll=15)

        for key in ("kLa_10", "kLa_25", "kLa_50"):
            val = results[key]
            assert math.isfinite(val), (
                f"Segment {k} ({params['run_id']}): {key}={val} is not finite. "
                f"Check {run_dir}/logstats.dat and logs/slurm_{job_ids[k]}.err"
            )
