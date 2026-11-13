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


def test_mpi_writes_its_own_record(tmp_path):
    src = _fake_src(tmp_path, {"BioReactor.c": "int main;"})
    seen = {}

    def run_case(case, binary, rd, params, restart, env):
        seen.update(case=case.name, ranks=case.ranks, level=params["fidelity"])
        _fake_run()(case, binary, rd, params, restart, env)

    rec = bb.run_mpi(src, tmp_path / "work", run_case=run_case,
                     build_fn=lambda src, work: work / "mpi", qoi_fn=lambda rd: {})
    on_disk = json.loads((tmp_path / "work" / "mpi4.json").read_text())
    assert on_disk == rec and rec["src_hash"] == bb.src_hash(src)
    assert rec["case"]["status"] == "ok" and rec["case"]["files"]
    assert seen == {"case": "mpi4", "ranks": 4, "level": 6}
    assert not (tmp_path / "work" / "bundle.json").exists()   # never touches the serial job's file


def test_restart_gets_the_parent_phase_sidecar(tmp_path):
    """The solver reads checkpoint.phase from its working directory (as submit_continue.stage_phase stages it)."""
    src = _fake_src(tmp_path, {"BioReactor.c": "int main;"})
    seen = {}

    def run_case(case, binary, rd, params, restart, env):
        if case.name == "kill":
            rd.mkdir(parents=True, exist_ok=True)
            (rd / "checkpoint.phase").write_text('{"t_phase_offset": 0}')
        if case.name == "kill_resume":
            seen["phase"] = (rd / "checkpoint.phase").read_text() if (rd / "checkpoint.phase").exists() else None
        _fake_run()(case, binary, rd, params, restart, env)

    bb.run_bundle(src, tmp_path / "work", only=["kill_resume"], run_case=run_case,
                  build_fn=lambda src, work, builds: {b: work / b for b in builds}, qoi_fn=lambda rd: {})
    assert seen["phase"] == '{"t_phase_offset": 0}'


def test_a_rerun_starts_from_an_empty_case_folder(tmp_path):
    """A file left by an earlier bundle in the same WORK_DIR must not be read as this run's output."""
    stale = tmp_path / "work" / "runs" / "cold" / "checkpoint.dump"
    stale.parent.mkdir(parents=True)
    stale.write_text("old")
    seen = {}

    def run_case(case, binary, rd, params, restart, env):
        if case.name == "cold":
            seen["stale"] = (rd / "checkpoint.dump").exists()
        _fake_run()(case, binary, rd, params, restart, env)

    src = _fake_src(tmp_path, {"BioReactor.c": "int main;"})
    bb.run_bundle(src, tmp_path / "work", only=["cold"], run_case=run_case,
                  build_fn=lambda src, work, builds: {b: work / b for b in builds}, qoi_fn=lambda rd: {})
    assert seen["stale"] is False
