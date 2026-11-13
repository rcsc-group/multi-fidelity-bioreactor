"""tools/bundle/bundle.py: runner contract, with a fake solver (no build, no run)."""
import json
import struct
from pathlib import Path

import pytest

from tools.bundle import bundle as bb
from tools.bundle import cases as bc


def _fake_src(d: Path, files: dict) -> Path:
    (d / "src").mkdir(parents=True)
    for name, text in files.items():
        (d / "src" / name).write_text(text)
    return d / "src"


def _fake_run(fail=()):
    def run_case(case, binary, rd, params, restart, env):
        if case.name in fail:
            raise RuntimeError(f"{case.name}: exit 1")
        rd.mkdir(parents=True, exist_ok=True)
        (rd / "normf.dat").write_text(case.name)
        (rd / "checkpoint.dump").write_bytes(struct.pack("<d", 1.7) + b"x" * 8)
        (rd / "replay_start.dump").write_bytes(struct.pack("<d", 1.0))
    return run_case


def _run(tmp_path, fail=()):
    src = _fake_src(tmp_path, {"BioReactor.c": "int main;", "a.h": "x"})
    return bb.run_bundle(src, tmp_path / "work", run_case=_fake_run(fail),
                         build_fn=lambda src, work, builds: {b: work / b for b in builds},
                         qoi_fn=lambda rd: {"tau_95": 1.0})


def test_failed_parent_marks_children(tmp_path):
    out = _run(tmp_path, fail={"half"})
    assert out["cases"]["half"]["status"] == "failed"
    assert "exit 1" in out["cases"]["half"]["error"]
    for child in ("cont", "caseB_warm"):
        assert out["cases"][child]["status"] == "not run"
        assert out["cases"][child]["error"] == "parent half failed"
    assert out["cases"]["cold"]["status"] == "ok"


def test_src_hash_ignores_file_order(tmp_path):
    a = _fake_src(tmp_path / "a", {"x.h": "1", "BioReactor.c": "2"})
    b = _fake_src(tmp_path / "b", {"BioReactor.c": "2", "x.h": "1"})
    c = _fake_src(tmp_path / "c", {"BioReactor.c": "2", "x.h": "changed"})
    assert bb.src_hash(a) == bb.src_hash(b) != bb.src_hash(c)


def test_bundle_json_schema(tmp_path):
    out = _run(tmp_path)
    on_disk = json.loads((tmp_path / "work" / "bundle.json").read_text())
    assert on_disk == out
    assert {"src_hash", "git_commit", "builds", "cases"} <= set(out)
    for name, c in out["cases"].items():
        assert {"status", "error", "wall_s", "dir", "files", "qoi"} <= set(c), name
    assert out["cases"]["cold"]["qoi"] == {"tau_95": 1.0}
    assert out["cases"]["mpi4"]["status"] == "not run"   # MPI runs in its own job
    assert "MPI" in out["cases"]["mpi4"]["error"]


def test_resume_takes_the_dump_time(tmp_path):
    dump = tmp_path / "checkpoint.dump"
    dump.write_bytes(struct.pack("<d", 1.7) + b"rest")
    p = bb.resume_params({"t_checkpoint": 9.9, "t_end": 9.9}, dump, t_end_abs=4.0)
    assert p["t_checkpoint"] == pytest.approx(1.7)
    assert p["t_end"] == pytest.approx(4.0 - 1.7)
