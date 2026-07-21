"""Tests for scripts/stage_segment.py — safe checkpoint staging for chain restarts.

Root cause this guards against: a chain-restart segment's checkpoint.dump can
exist on disk (e.g. copied in as an unrun restart input, or left over from a
job that was killed seconds after starting) without the segment ever having
actually finished. Only a well-formed results.json proves a segment ran to
completion, since it's written exclusively by postprocess.py at the end of a
real run. stage_segment.py must refuse to stage a "next" segment from a
"prev" segment's checkpoint unless prev's results.json says so.
"""
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parents[1]))
from scripts.stage_segment import SegmentNotCompleteError, stage_segment


def _make_run(runs_root: Path, run_id: str, *, with_results: bool, with_checkpoint: bool = True,
              with_params: bool = True) -> Path:
    run_dir = runs_root / run_id
    run_dir.mkdir(parents=True)
    if with_params:
        (run_dir / "params.json").write_text(json.dumps({"run_id": run_id}))
    if with_checkpoint:
        (run_dir / "checkpoint.dump").write_bytes(b"fake-dump-data")
    if with_results:
        (run_dir / "results.json").write_text(json.dumps({"tau_100_max": 0.123}))
    return run_dir


def test_refuses_to_stage_when_prev_has_no_results_json(tmp_path):
    """This is the exact bug: prev's checkpoint.dump exists (stale/copied-in),
    but prev never actually finished (no results.json). Must refuse."""
    runs_root = tmp_path / "runs"
    scratch_root = tmp_path / "scratch"
    _make_run(runs_root, "prev_rid", with_results=False)
    _make_run(runs_root, "next_rid", with_results=False)

    with pytest.raises(SegmentNotCompleteError):
        stage_segment("prev_rid", "next_rid", runs_root=runs_root, scratch_root=scratch_root)

    assert not (scratch_root / "next_rid").exists()


def test_stages_successfully_when_prev_is_genuinely_complete(tmp_path):
    runs_root = tmp_path / "runs"
    scratch_root = tmp_path / "scratch"
    _make_run(runs_root, "prev_rid", with_results=True)
    _make_run(runs_root, "next_rid", with_results=False)

    stage_segment("prev_rid", "next_rid", runs_root=runs_root, scratch_root=scratch_root)

    staged_ckpt = scratch_root / "next_rid" / "checkpoint.dump"
    staged_params = scratch_root / "next_rid" / "params.json"
    assert staged_ckpt.exists()
    assert staged_ckpt.read_bytes() == (runs_root / "prev_rid" / "checkpoint.dump").read_bytes()
    assert staged_params.exists()


def test_refuses_when_prev_results_json_is_malformed(tmp_path):
    runs_root = tmp_path / "runs"
    scratch_root = tmp_path / "scratch"
    prev_dir = _make_run(runs_root, "prev_rid", with_results=False)
    (prev_dir / "results.json").write_text("{not valid json")
    _make_run(runs_root, "next_rid", with_results=False)

    with pytest.raises(SegmentNotCompleteError):
        stage_segment("prev_rid", "next_rid", runs_root=runs_root, scratch_root=scratch_root)


def test_refuses_when_next_params_json_missing(tmp_path):
    """Staging needs the next segment's own params.json (chain-generated) to exist."""
    runs_root = tmp_path / "runs"
    scratch_root = tmp_path / "scratch"
    _make_run(runs_root, "prev_rid", with_results=True)
    _make_run(runs_root, "next_rid", with_results=False, with_params=False, with_checkpoint=False)

    with pytest.raises(FileNotFoundError):
        stage_segment("prev_rid", "next_rid", runs_root=runs_root, scratch_root=scratch_root)


def test_refuses_when_prev_checkpoint_missing(tmp_path):
    runs_root = tmp_path / "runs"
    scratch_root = tmp_path / "scratch"
    _make_run(runs_root, "prev_rid", with_results=True, with_checkpoint=False)
    _make_run(runs_root, "next_rid", with_results=False)

    with pytest.raises(FileNotFoundError):
        stage_segment("prev_rid", "next_rid", runs_root=runs_root, scratch_root=scratch_root)
