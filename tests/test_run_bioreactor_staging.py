"""tests/conftest.run_bioreactor stages a restart the way production does (scripts/submit_continue.stage_phase):
the dump's label, checkpoint.phase, goes next to the dump in the run folder. Without it the solver refuses the
restart (src/forcing_state.h restart_source), so every medium test that restarts would fail in CI."""
import json

import tests.conftest as cf


def _fake_binary(tmp_path):
    """A stand-in solver: lists the files of its working directory into seen.txt."""
    b = tmp_path / "fake_solver.sh"
    b.write_text("#!/bin/sh\nls > seen.txt\n")
    b.chmod(0o755)
    return b


def test_restart_stages_the_dumps_label(tmp_path, monkeypatch):
    monkeypatch.setattr(cf, "BINARY", _fake_binary(tmp_path))
    parent = tmp_path / "parent"
    parent.mkdir()
    (parent / "checkpoint.dump").write_bytes(b"\0" * 16)
    (parent / "checkpoint.phase").write_text(json.dumps({"t": 1.0}))
    rd = cf.run_bioreactor({"run_id": "child"}, tmp_path / "runs", restart_from=parent / "checkpoint.dump")
    seen = (rd / "seen.txt").read_text().split()
    assert "checkpoint.phase" in seen and "restart.dump" in seen
    assert (rd / "checkpoint.phase").read_text() == (parent / "checkpoint.phase").read_text()


def test_restart_without_a_label_stages_nothing(tmp_path, monkeypatch):
    """A dump with no label (the refusal path under test) must stay without one."""
    monkeypatch.setattr(cf, "BINARY", _fake_binary(tmp_path))
    parent = tmp_path / "parent"
    parent.mkdir()
    (parent / "checkpoint.dump").write_bytes(b"\0" * 16)
    rd = cf.run_bioreactor({"run_id": "child"}, tmp_path / "runs", restart_from=parent / "checkpoint.dump")
    assert "checkpoint.phase" not in (rd / "seen.txt").read_text().split()
