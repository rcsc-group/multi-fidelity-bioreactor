"""scripts/paths.py: one place for machine-specific folders (environment variable, OSCAR default)."""
import importlib
from pathlib import Path

from scripts import paths


def test_repo_is_this_checkout():
    assert paths.REPO == Path(__file__).resolve().parents[1]


def test_oscar_defaults_without_environment(monkeypatch):
    for v in ("MFB_SCRATCH", "MFB_DATA", "BASILISK", "MFB_LEGACY_2D"):
        monkeypatch.delenv(v, raising=False)
    p = importlib.reload(paths)
    assert str(p.SCRATCH) == "/oscar/scratch/eaguerov"
    assert str(p.DATA) == "/oscar/data/dharri15/eaguerov"
    assert str(p.BASILISK) == "/oscar/data/dharri15/eaguerov/basilisk/src"


def test_environment_overrides_the_default(monkeypatch, tmp_path):
    monkeypatch.setenv("MFB_SCRATCH", str(tmp_path))
    p = importlib.reload(paths)
    assert p.SCRATCH == tmp_path
    assert f"{p.SCRATCH}/mpi_runs" == f"{tmp_path}/mpi_runs"
    monkeypatch.delenv("MFB_SCRATCH")
    importlib.reload(paths)


def test_expand_replaces_known_variables_only(monkeypatch, tmp_path):
    monkeypatch.setenv("MFB_SCRATCH", str(tmp_path))
    p = importlib.reload(paths)
    assert p.expand("$MFB_SCRATCH/BioReactor-mpi") == f"{tmp_path}/BioReactor-mpi"
    assert p.expand("${MFB_SCRATCH}/x") == f"{tmp_path}/x"
    assert p.expand("$HOME/x") == "$HOME/x"
    monkeypatch.delenv("MFB_SCRATCH")
    importlib.reload(paths)
