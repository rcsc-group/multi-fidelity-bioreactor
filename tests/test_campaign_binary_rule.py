"""Campaign 2D: a warm source must come from the same solver source as the restarting run; other sources are refused.
The 2026-11-04 rule compared binary sha256; since 2026-11-05 it compares the build sidecar <binary>.build.json
(commit + flags; see test_campaign_source_rule.py). Synthetic binaries, each with a sidecar, and run folders only.
"""
import hashlib
import json
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "campaign2d"))

import campaign as c  # noqa: E402


def _bin(path: Path, commit: str) -> Path:
    """A synthetic binary plus its build sidecar."""
    path.write_bytes(commit.encode())
    Path(str(path) + ".build.json").write_text(json.dumps({"commit": commit, "flags": ["-D_MPI=1"]}))
    return path


def _src(root: Path, name: str, binary: Path, rpm=20.0, theta=4.0):
    d = root / name
    d.mkdir(parents=True)
    p = {"fidelity": 8, "omega_b": rpm * 2 * math.pi / 60.0, "theta_max": [theta, 0.0, 0.0], "_binary": str(binary)}
    (d / "params.json").write_text(json.dumps(p))
    (d / "results.json").write_text("{}")
    (d / "checkpoint.dump").write_bytes(b"x")
    return name


def test_same_binary_source_is_accepted(tmp_path):
    target = _bin(tmp_path / "bin_2ef5b07", "2ef5b07")
    name = _src(tmp_path / "runs", "src_same", target)
    smap = {(20.0, 4.0): name}
    assert c.require_warm_sources([[20.0, 4.0]], smap, target, runs_root=tmp_path / "runs") == [name]


def test_cross_binary_source_is_refused(tmp_path):
    target = _bin(tmp_path / "bin_2ef5b07", "2ef5b07")
    other = _bin(tmp_path / "bin_f1c11e0", "f1c11e0")
    name = _src(tmp_path / "runs", "src_other", other)
    with pytest.raises(ValueError, match="binary"):
        c.require_warm_sources([[20.0, 4.0]], {(20.0, 4.0): name}, target, runs_root=tmp_path / "runs")


def test_source_with_missing_binary_is_refused(tmp_path):
    target = _bin(tmp_path / "bin_2ef5b07", "2ef5b07")
    gone = tmp_path / "bin_purged"          # not on disk: hash unknown, so not the same binary
    name = _src(tmp_path / "runs", "src_gone", gone)
    with pytest.raises(ValueError, match="binary"):
        c.require_warm_sources([[20.0, 4.0]], {(20.0, 4.0): name}, target, runs_root=tmp_path / "runs")


def test_binary_hash_is_sha256(tmp_path):
    f = tmp_path / "b"
    f.write_bytes(b"abc")
    assert c.binary_sha256(f) == hashlib.sha256(b"abc").hexdigest()


def test_same_binary_map_keeps_only_matching_sources(tmp_path):
    target = _bin(tmp_path / "bin_2ef5b07", "2ef5b07")
    other = _bin(tmp_path / "bin_f1c11e0", "f1c11e0")
    _src(tmp_path / "runs", "good", target, rpm=20.0, theta=4.0)
    _src(tmp_path / "runs", "bad", other, rpm=25.0, theta=4.0)
    m = c.same_binary_map(target, runs_root=tmp_path / "runs")
    assert m == {(20.0, 4.0): "good"}
