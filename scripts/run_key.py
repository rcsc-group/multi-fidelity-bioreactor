"""Run keys: find a finished run instead of computing it again.

A run key is the sha256 of three things:

1. the physics in params.json: every key except ``run_id`` and the bookkeeping keys that start with ``_``.
   Floats are rounded to 12 significant digits, and 5.0 is the same as 5;
2. the binary: the sha256 of the binary file (``binary_identity``);
3. the parent: the key of the run whose checkpoint.dump a restart reads (None for a cold start).

Two runs with the same key give the same results. ``simulate.submit_slurm`` writes the key into params.json
(``_run_key``) and into the index, and refuses a run whose key is already in the index (``DuplicateRun``).

The index is ``<runs_root>/_index.jsonl``: one JSON line per run, the last line of a run_id wins. Submissions
append to it; ``python -m scripts.run_key index [RUNS_ROOT]`` rebuilds it from the run folders (run it as a
SLURM job: it reads every params.json and hashes every binary).

A key is usable for duplicate checks (``dedupe``) only if the binary is known and, for a restart, the parent's key
is usable. A binary is unknown if its file is gone, or if it changed after the run (rebuilt in place under the same
name), because then nobody knows which code made the run.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # repo root
from scripts import paths  # noqa: E402

INDEX_NAME = "_index.jsonl"
SOURCE_KEYS = ("_parent_run", "_source_run")   # where older restarts recorded the run they started from
_SIG = 12                                       # significant digits kept in floats


def _canon(v):
    if isinstance(v, bool) or v is None or isinstance(v, str):
        return v
    if isinstance(v, (int, float)):
        f = float(f"{float(v):.{_SIG}g}")
        return int(f) if f.is_integer() else f
    if isinstance(v, dict):
        return {k: _canon(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_canon(x) for x in v]
    return str(v)


def physics(params: dict) -> dict:
    """The params that set the results: all but run_id and the ``_`` bookkeeping keys, in canonical form."""
    return _canon({k: v for k, v in params.items() if k != "run_id" and not k.startswith("_")})


def run_key(params: dict, binary_id: str | None, parent_key: str | None) -> str:
    """sha256 (hex) of the physics, the binary identity and the parent run's key."""
    blob = json.dumps({"physics": physics(params), "binary": binary_id, "parent": parent_key},
                      sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode()).hexdigest()


_SHA_CACHE: dict = {}


def binary_identity(binary, not_after: float | None = None) -> str | None:
    """``sha256:<hex>`` of the binary file. None (unknown) if there is no file, or if the file changed after
    ``not_after`` (a time.time() value: the run used an older build with the same name)."""
    if not binary or not Path(binary).is_file():
        return None
    st = Path(binary).stat()
    if not_after is not None and st.st_mtime > not_after:
        return None
    ck = (str(binary), st.st_mtime, st.st_size)
    if ck not in _SHA_CACHE:
        _SHA_CACHE[ck] = "sha256:" + hashlib.sha256(Path(binary).read_bytes()).hexdigest()
    return _SHA_CACHE[ck]


def default_binary(template, project_root) -> str:
    """The binary a SLURM template runs when params.json has no ``_binary`` (see the BINARY lines there)."""
    if Path(template).name == "slurm_mpi_template.sh":
        return f"{paths.SCRATCH}/BioReactor-mpi-video"
    return str(Path(project_root) / "build" / "BioReactor-video")


# -- the index ---------------------------------------------------------------------------------------------------

def index_path(runs_root) -> Path:
    return Path(runs_root) / INDEX_NAME


def read_index(runs_root) -> list[dict]:
    """The index entries, one per run_id (the last line of a run_id wins). Empty if there is no index."""
    p = index_path(runs_root)
    if not p.exists():
        return []
    out = {}
    for line in p.read_text().splitlines():
        if line.strip():
            e = json.loads(line)
            out[e["run_id"]] = e
    return list(out.values())


def append_index(runs_root, entry: dict) -> None:
    """Add one entry (one line, one write call, so concurrent submitters do not mix lines)."""
    p = index_path(runs_root)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "a") as f:
        f.write(json.dumps(entry, sort_keys=True) + "\n")


def entry(run_id, params, binary_id, parent_run, parent_entry) -> dict:
    """One index entry. ``parent_entry`` is the parent's index entry (None if unknown or a cold start)."""
    restart = float(params.get("t_checkpoint", 0) or 0) > 0
    parent_key = parent_entry["key"] if parent_entry else None
    dedupe = binary_id is not None and (not restart or bool(parent_entry and parent_entry["dedupe"]))
    return {"run_id": run_id, "key": run_key(params, binary_id, parent_key), "binary": binary_id,
            "parent_run": parent_run, "parent_key": parent_key, "dedupe": dedupe, "time": time.time()}


def job_alive(job: str) -> bool:
    """True if SLURM still has the job (pending or running)."""
    import subprocess
    r = subprocess.run(["squeue", "-h", "-j", str(job), "-o", "%T"], capture_output=True, text=True)
    return r.returncode == 0 and bool(r.stdout.strip())


def find_duplicate(runs_root, e: dict) -> dict | None:
    """An indexed run with the same usable key as entry ``e`` and another run_id that is finished (results.json)
    or still in SLURM, or None. A failed run (no results, job gone) does not count: the point may be retried."""
    if not e["dedupe"]:
        return None
    for o in read_index(runs_root):
        if not (o["dedupe"] and o["key"] == e["key"] and o["run_id"] != e["run_id"]):
            continue
        if (Path(runs_root) / o["run_id"] / "results.json").exists() or (o.get("job") and job_alive(o["job"])):
            return o
    return None


def build_index(runs_root) -> list[dict]:
    """Key every run folder under runs_root that has a params.json and write a new index. Returns the entries.

    Old runs: the binary is identified now (unknown if it changed after the run's params.json was written); the
    parent of a restart comes from _parent_run or _source_run, else it is unknown."""
    root = Path(runs_root)
    params = {}
    for d in sorted(p for p in root.iterdir() if p.is_dir()):
        f = d / "params.json"
        if f.is_file():
            try:
                params[d.name] = (json.loads(f.read_text()), f.stat().st_mtime)
            except (json.JSONDecodeError, OSError):
                continue
    done: dict = {}

    def key_of(rid, stack=()):
        if rid in done:
            return done[rid]
        if rid not in params or rid in stack:
            return None
        p, mtime = params[rid]
        parent = next((p[k] for k in SOURCE_KEYS if p.get(k)), None)
        restart = float(p.get("t_checkpoint", 0) or 0) > 0
        parent_entry = key_of(parent, stack + (rid,)) if (restart and parent) else None
        done[rid] = entry(rid, p, binary_identity(p.get("_binary"), not_after=mtime),
                          parent if restart else None, parent_entry)
        return done[rid]

    for rid in params:
        key_of(rid)
    entries = [done[r] for r in params]
    tmp = index_path(root).with_suffix(".tmp")
    tmp.write_text("".join(json.dumps(e, sort_keys=True) + "\n" for e in entries))
    os.replace(tmp, index_path(root))
    return entries


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if not argv or argv[0] != "index":
        print(__doc__)
        return 2
    root = Path(argv[1]) if len(argv) > 1 else Path(__file__).resolve().parents[1] / "runs"
    entries = build_index(root)
    n_dedupe = sum(e["dedupe"] for e in entries)
    print(f"{len(entries)} runs indexed in {index_path(root)}; {n_dedupe} with a usable key")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
