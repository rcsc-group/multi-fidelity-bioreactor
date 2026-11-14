"""Which solver source made a run, and whether two runs form a clean pair for a feature's effect.

A run names its binary in params.json (`_binary`). The binary's sidecar `<binary>.build.json` gives the source
commit and the compile flags. Without a sidecar, a binary whose sha256 equals that of a binary with a sidecar
has the same source ("sha-match"). Anything else is "unknown": the ledger counts such runs but never uses them
for a number, because their source cannot be named.
"""
import hashlib
import json
from pathlib import Path

from scripts.campaign2d.campaign import build_info

BIN_DIR = Path("/oscar/data/dharri15/eaguerov/bin")


def _sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def known_builds(bin_dir: Path = BIN_DIR) -> dict:
    """sha256 -> sidecar dict, for every binary in bin_dir that has a `.build.json` sidecar."""
    out = {}
    for side in sorted(Path(bin_dir).glob("*.build.json")):
        info = json.loads(side.read_text())
        b = Path(str(side)[: -len(".build.json")])
        sha = info.get("sha256") or (_sha256(b) if b.is_file() else None)
        if sha:
            out[sha] = info
    return out


def _on_disk(binary: str | None, bin_dir: Path) -> Path | None:
    """The binary, or its persistent copy in bin_dir with the same file name (scratch purges unused files)."""
    if not binary:
        return None
    if Path(binary).is_file():
        return Path(binary)
    alt = Path(bin_dir) / Path(binary).name
    return alt if alt.is_file() else None


def provenance(run_dir: Path, known: dict, bin_dir: Path = BIN_DIR) -> dict:
    """{"commit", "flags", "status": "sidecar"|"sha-match"|"unknown", "sha"} of the binary that made run_dir
    ("sha" is None when the binary is no longer on disk)."""
    binary = json.loads((Path(run_dir) / "params.json").read_text()).get("_binary")
    found = _on_disk(binary, bin_dir)
    sha = _sha256(found) if found else None
    info = build_info(binary)
    if info is not None:
        return {"commit": info.get("commit"), "flags": list(info.get("flags", [])), "status": "sidecar", "sha": sha}
    info = known.get(sha) if sha else None
    if info is not None:
        return {"commit": info.get("commit"), "flags": list(info.get("flags", [])), "status": "sha-match", "sha": sha}
    return {"commit": None, "flags": [], "status": "unknown", "sha": sha}


def pair_status(a: dict, b: dict, mechanism_flag: str | None) -> str:
    """"same-source" if both runs come from one known commit and their flags differ at most by the mechanism's
    own compile flag, or both ran byte-identical binaries (same sha256, commit unknown); else "confounded" (the
    difference mixes the feature with a source change)."""
    if a.get("sha") and a.get("sha") == b.get("sha") and not mechanism_flag:
        return "same-source"
    if a["status"] == "unknown" or b["status"] == "unknown" or a["commit"] != b["commit"]:
        return "confounded"
    diff = set(a["flags"]) ^ set(b["flags"])
    return "same-source" if diff <= ({mechanism_flag} if mechanism_flag else set()) else "confounded"


RUNS = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/runs")
PAIRS = Path(__file__).with_name("pairs.yaml")


def load_pairs(path: Path = PAIRS, runs: Path = RUNS) -> dict:
    """pairs.yaml with every run name resolved to a folder (absolute names kept, others under `runs`)."""
    import yaml
    out = yaml.safe_load(Path(path).read_text())
    for mech in out.values():
        for p in mech["pairs"]:
            for k in ("with", "without"):
                p[k] = Path(p[k]) if Path(p[k]).is_absolute() else Path(runs) / p[k]
    return out
