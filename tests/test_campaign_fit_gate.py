"""Campaign 2D (2026-11-04): FAILED-GATE rows leave the fit; every warm submission carries its provenance.

Synthetic data only. A FAILED-GATE run is dropped; PASS, PENDING and cold runs (no gate entry) stay.
A warm submission without _source_run or _err_file is refused.
"""
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "campaign2d"))

import campaign as c  # noqa: E402


def _db():
    return pd.DataFrame({"level": [8, 8, 9, 10], "rpm": [20.0, 25.0, 20.0, 30.0], "theta": [3.0, 3.0, 3.0, 4.0],
                         "y": [1.0, 2.0, 3.0, 4.0], "run": ["a", "b", "c", "d"]})


def test_failed_gate_rows_are_dropped():
    out = c.exclude_failed_gate(_db(), {"b": "FAILED-GATE", "c": "PASS", "d": "PENDING"})
    assert list(out.run) == ["a", "c", "d"]


def test_cold_runs_without_gate_entry_are_kept():
    out = c.exclude_failed_gate(_db(), {"b": "FAILED-GATE"})
    assert "a" in set(out.run)
    assert len(out) == 3


def test_warm_submission_refused_without_source_run():
    with pytest.raises(ValueError, match="_source_run"):
        c.warm_provenance({"run_id": "x", "t_checkpoint": 5.0}, source_run=None, err_file="/x/slurm.err")


def test_warm_submission_refused_without_err_file():
    with pytest.raises(ValueError, match="_err_file"):
        c.warm_provenance({"run_id": "x", "t_checkpoint": 5.0}, source_run="src", err_file=None)


def test_cold_submission_refused_as_warm():
    with pytest.raises(ValueError, match="t_checkpoint"):
        c.warm_provenance({"run_id": "x", "t_checkpoint": 0.0}, source_run="src", err_file="/x/slurm.err")


def test_warm_provenance_sets_both_keys():
    q = c.warm_provenance({"run_id": "x", "t_checkpoint": 5.0}, source_run="src", err_file="/x/slurm.err")
    assert q["_source_run"] == "src" and q["_err_file"] == "/x/slurm.err" and q["run_id"] == "x"
