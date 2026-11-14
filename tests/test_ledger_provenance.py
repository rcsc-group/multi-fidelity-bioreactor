"""tools/ledger/provenance.py: which source made a run, and whether two runs form a clean (same-source) pair."""
import hashlib
import json
from pathlib import Path

from tools.ledger import provenance as pv


def _binary(tmp: Path, name: str, body: bytes, sidecar: dict | None = None) -> Path:
    b = tmp / name
    b.write_bytes(body)
    if sidecar is not None:
        Path(str(b) + ".build.json").write_text(json.dumps(sidecar))
    return b


def _run(tmp: Path, name: str, binary: Path) -> Path:
    d = tmp / name
    d.mkdir()
    (d / "params.json").write_text(json.dumps({"_binary": str(binary)}))
    return d


def test_sidecar_wins(tmp_path):
    b = _binary(tmp_path, "bin_a", b"A", {"commit": "abc1234", "flags": ["-D_MPI=1"]})
    known = {hashlib.sha256(b"A").hexdigest(): {"commit": "other99", "flags": []}}
    p = pv.provenance(_run(tmp_path, "r", b), known)
    assert {k: p[k] for k in ("commit", "flags", "status")} == {"commit": "abc1234", "flags": ["-D_MPI=1"], "status": "sidecar"}


def test_sha_match_without_sidecar(tmp_path):
    b = _binary(tmp_path, "copy_of_a", b"A")
    known = {hashlib.sha256(b"A").hexdigest(): {"commit": "abc1234", "flags": ["-D_MPI=1"]}}
    p = pv.provenance(_run(tmp_path, "r", b), known)
    assert {k: p[k] for k in ("commit", "flags", "status")} == {"commit": "abc1234", "flags": ["-D_MPI=1"], "status": "sha-match"}


def test_unknown_binary_is_counted_not_used(tmp_path):
    b = _binary(tmp_path, "mystery", b"Z")
    p = pv.provenance(_run(tmp_path, "r", b), {})
    assert p["status"] == "unknown" and p["commit"] is None
    missing = pv.provenance(_run(tmp_path, "r2", tmp_path / "gone"), {}, bin_dir=tmp_path / "nobin")
    assert missing["status"] == "unknown" and missing["sha"] is None
    assert pv.pair_status(missing, missing, None) == "confounded"   # no commit, no binary: never a clean pair


def test_cross_source_pair_is_confounded():
    a = {"commit": "abc1234", "flags": ["-D_MPI=1"], "status": "sidecar"}
    b = {"commit": "def5678", "flags": ["-D_MPI=1"], "status": "sidecar"}
    assert pv.pair_status(a, b, None) == "confounded"
    c = {"commit": "abc1234", "flags": ["-D_MPI=1", "-DEXTRA_TRACERS=1"], "status": "sidecar"}
    assert pv.pair_status(a, c, None) == "confounded"   # a flag that is not the mechanism's


def test_cross_level_flag_only_is_same_source():
    a = {"commit": "abc1234", "flags": ["-D_MPI=1"], "status": "sidecar"}
    b = {"commit": "abc1234", "flags": ["-D_MPI=1", "-DCROSS_LEVEL_WARMSTART=1"], "status": "sha-match"}
    assert pv.pair_status(a, b, "-DCROSS_LEVEL_WARMSTART=1") == "same-source"
    assert pv.pair_status(a, b, None) == "confounded"
    assert pv.pair_status(a, dict(a), None) == "same-source"


def test_pairs_file_names_existing_runs():
    """Every run in pairs.yaml is on disk with a params.json."""
    pairs = pv.load_pairs()
    assert set(pairs) == {"restart_label", "cross_level_warmstart", "chaining", "flux_replay"}
    for mech, m in pairs.items():
        assert m["qoi"], mech
        for p in m["pairs"]:
            for k in ("with", "without"):
                assert (p[k] / "params.json").is_file(), f"{mech}: {p[k]} has no params.json"


def test_same_unknown_binary_is_same_source(tmp_path):
    """Two runs on byte-identical binaries share a source even when no sidecar names the commit."""
    b1 = _binary(tmp_path, "x", b"SAME")
    b2 = _binary(tmp_path, "x_copy", b"SAME")
    b3 = _binary(tmp_path, "y", b"OTHER")
    p1, p2, p3 = (pv.provenance(_run(tmp_path, f"r{i}", b), {}) for i, b in enumerate((b1, b2, b3)))
    assert p1["status"] == "unknown"
    assert pv.pair_status(p1, p2, None) == "same-source"
    assert pv.pair_status(p1, p3, None) == "confounded"


def test_purged_binary_is_found_by_name_in_bin_dir(tmp_path):
    """A binary purged from scratch resolves to the persistent copy with the same file name."""
    (tmp_path / "bin").mkdir()
    kept = _binary(tmp_path / "bin", "BioReactor-mpi-xlevel-chain", b"X")
    other = _binary(tmp_path, "x_other_path", b"X")
    gone = tmp_path / "scratch" / "BioReactor-mpi-xlevel-chain"
    a = pv.provenance(_run(tmp_path, "r1", gone), {}, bin_dir=tmp_path / "bin")
    b = pv.provenance(_run(tmp_path, "r2", other), {}, bin_dir=tmp_path / "bin")
    assert a["sha"] == b["sha"] is not None
    assert pv.pair_status(a, b, None) == "same-source"
    assert kept.is_file()
