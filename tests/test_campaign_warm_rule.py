"""Campaign 2D (lead rule, 2026-11-04): every L9 and L10 point warm-starts from a COMPLETED L8 run at the same
(rpm, theta). Cold points are refused. Synthetic run folders only.
"""
import json
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "campaign2d"))

import campaign as c  # noqa: E402


def _run(root: Path, name: str, rpm: float, theta: float, fidelity: int = 8, results=True, checkpoint=True):
    d = root / name
    d.mkdir(parents=True)
    p = {"fidelity": fidelity, "omega_b": rpm * 2 * math.pi / 60.0, "theta_max": [theta, 0.0, 0.0], "_binary": __file__}
    (d / "params.json").write_text(json.dumps(p))
    if results:
        (d / "results.json").write_text("{}")
    if checkpoint:
        (d / "checkpoint.dump").write_bytes(b"x")
    return d


def test_only_completed_L8_runs_are_sources(tmp_path):
    _run(tmp_path, "done_l8", 20.0, 4.0)                       # completed L8: a source
    _run(tmp_path, "no_ck_l8", 25.0, 4.0, checkpoint=False)    # no checkpoint: not a source
    _run(tmp_path, "unfinished", 30.0, 4.0, results=False)     # not finished: not a source
    _run(tmp_path, "l9_run", 22.5, 5.0, fidelity=9)            # L9 run: not an L8 source
    smap = c.l8_source_map(runs_root=tmp_path)
    assert smap == {(20.0, 4.0): "done_l8"}


def test_cold_point_is_refused():
    smap = {(20.0, 4.0): "done_l8"}
    with pytest.raises(ValueError, match="cold"):
        c.require_warm_sources([[25.0, 4.0]], smap, __file__)


def test_warm_point_gets_its_source_run(tmp_path):
    # sources and target share one binary with a build sidecar (same-source rule, test_campaign_source_rule.py)
    binary = tmp_path / "bin"
    binary.write_bytes(b"x")
    (tmp_path / "bin.build.json").write_text(json.dumps({"commit": "2ef5b07", "flags": ["-D_MPI=1"]}))
    runs = tmp_path / "runs"
    for name, rpm, th in (("done_l8", 20.0, 4.0), ("other_l8", 30.0, 6.0)):
        d = _run(runs, name, rpm, th)
        p = json.loads((d / "params.json").read_text())
        p["_binary"] = str(binary)
        (d / "params.json").write_text(json.dumps(p))
    smap = {(20.0, 4.0): "done_l8", (30.0, 6.0): "other_l8"}
    got = c.require_warm_sources([[20.0, 4.0], [30.0, 6.0]], smap, binary, runs_root=runs)
    assert got == ["done_l8", "other_l8"]


def test_measured_warm_cost_is_the_recorded_median():
    assert c.WARM_MEDIAN_CORE_H[9] == pytest.approx(97.7)
    assert c.NEXT_COST_L9 == c.WARM_MEDIAN_CORE_H[9]
