"""Campaign 2D (lead rule, 2026-11-05): a warm source is valid when its binary has the SAME SOURCE COMMIT as the
restarting binary, and the flags differ only by ALLOWED_EXTRA_FLAGS (-DCROSS_LEVEL_WARMSTART=1). Build sidecar
<binary>.build.json holds {commit, flags}. Grid guard: a warm restart must have total cells == 4**level (uniform).
Synthetic files only.
"""
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "campaign2d"))

import campaign as c  # noqa: E402


def _bin(root: Path, name: str, commit: str, flags):
    b = root / name
    b.write_bytes(name.encode())
    (root / (name + ".build.json")).write_text(json.dumps({"commit": commit, "flags": flags}))
    return b


def _src(root: Path, name: str, binary: Path, rpm=20.0, theta=4.0):
    d = root / "runs" / name
    d.mkdir(parents=True)
    p = {"fidelity": 8, "omega_b": rpm * 2 * 3.141592653589793 / 60.0, "theta_max": [theta, 0.0, 0.0], "_binary": str(binary)}
    (d / "params.json").write_text(json.dumps(p))
    (d / "results.json").write_text("{}")
    (d / "checkpoint.dump").write_bytes(b"x")
    return name


def test_same_commit_cross_flag_is_same_source(tmp_path):
    lean = _bin(tmp_path, "lean", "2ef5b07", ["-D_MPI=1", "-DEXTRA_TRACERS=0"])
    xl = _bin(tmp_path, "xl", "2ef5b07", ["-D_MPI=1", "-DEXTRA_TRACERS=0", "-DCROSS_LEVEL_WARMSTART=1"])
    assert c.same_source(lean, xl)
    assert c.same_source(xl, lean)


def test_different_commit_is_not_same_source(tmp_path):
    a = _bin(tmp_path, "a", "2ef5b07", ["-D_MPI=1"])
    b = _bin(tmp_path, "b", "f1c11e0", ["-D_MPI=1"])
    assert not c.same_source(a, b)


def test_unknown_extra_flag_is_not_same_source(tmp_path):
    a = _bin(tmp_path, "a", "2ef5b07", ["-D_MPI=1"])
    b = _bin(tmp_path, "b", "2ef5b07", ["-D_MPI=1", "-DDIAGNOSTICS=1"])
    assert not c.same_source(a, b)


def test_missing_sidecar_is_not_same_source(tmp_path):
    a = _bin(tmp_path, "a", "2ef5b07", ["-D_MPI=1"])
    b = tmp_path / "nosidecar"
    b.write_bytes(b"x")
    assert not c.same_source(a, b)


def test_require_refuses_other_commit_accepts_cross_flag(tmp_path):
    target = _bin(tmp_path, "target", "2ef5b07", ["-D_MPI=1", "-DCROSS_LEVEL_WARMSTART=1"])
    good = _bin(tmp_path, "good", "2ef5b07", ["-D_MPI=1"])
    bad = _bin(tmp_path, "bad", "f1c11e0", ["-D_MPI=1"])
    runs = tmp_path / "runs"
    _src(tmp_path, "src_good", good, rpm=20.0, theta=4.0)
    _src(tmp_path, "src_bad", bad, rpm=25.0, theta=4.0)
    assert c.require_warm_sources([[20.0, 4.0]], {(20.0, 4.0): "src_good"}, target, runs_root=runs) == ["src_good"]
    with pytest.raises(ValueError, match="source"):
        c.require_warm_sources([[25.0, 4.0]], {(25.0, 4.0): "src_bad"}, target, runs_root=runs)


def test_same_source_map_keeps_only_same_commit(tmp_path):
    target = _bin(tmp_path, "target", "2ef5b07", ["-D_MPI=1", "-DCROSS_LEVEL_WARMSTART=1"])
    good = _bin(tmp_path, "good", "2ef5b07", ["-D_MPI=1"])
    bad = _bin(tmp_path, "bad", "f1c11e0", ["-D_MPI=1"])
    _src(tmp_path, "g", good, rpm=20.0, theta=4.0)
    _src(tmp_path, "b", bad, rpm=25.0, theta=4.0)
    m = c.same_source_map(target, runs_root=tmp_path / "runs")
    assert m == {(20.0, 4.0): "g"}


def test_grid_guard_uniform_total_cells():
    assert c.grid_ok(total_cells=65536, level=8)          # L8: 256^2
    assert c.grid_ok(total_cells=262144, level=9)         # L9: 512^2
    assert not c.grid_ok(total_cells=32 * 2048, level=9)  # 2048 per rank x 32: the L8 grid, refused
    assert not c.grid_ok(total_cells=16 * 4096, level=9)


def test_parse_first_cells_line():
    line = "i: 1 t: 0.01 dt: 1e-4 #Cells: 2048 Wall clock time (s): 1.2 CPU time (s): 1.1"
    assert c.first_cells_per_rank(line + "\n") == 2048
    assert c.first_cells_per_rank("no cells here\n") is None
