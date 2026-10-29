"""Divergence of the face flux that a record feeds the scalars: the true flux (lfrec_l8_*) against the Kim-style rebuild (kimrec_l8_*).
Submits two 1-cycle replays (replay_diag = 1) from the existing records; each writes replay_diag.txt (src/BioReactor.c, replay_diag):
per step `iter t dfar dall umax`, with d = |sum over the 4 faces of the cell of the outward flux difference| (velocity units, = Delta*div),
dfar over cells with cs = 1 and four open faces, dall over all cells. The metric is d / umax.
Usage: uv run python scripts/diag_div_record.py submit | report
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import replay_gate as rg  # noqa: E402
from scripts import replay_lf as rl  # noqa: E402
from scripts import postprocess as pp  # noqa: E402

RPM, THETA = 32.5, 7.0
CASES = {  # tag -> (record run id, record dir)
    "true": (f"lfrec_l8_rpm{RPM:g}_th{THETA:g}", Path("/oscar/data/dharri15/eaguerov/replay_records/lf") / f"lf_l8_rpm{RPM:g}_th{THETA:g}"),
    "kim": (f"kimrec_l8_rpm{RPM:g}_th{THETA:g}", Path("/oscar/data/dharri15/eaguerov/replay_records/kim") / f"kim_l8_rpm{RPM:g}_th{THETA:g}"),
}


def diag_id(tag: str) -> str:
    return f"divdiag_{tag}_l8_rpm{RPM:g}_th{THETA:g}"


def submit() -> None:
    src = rl.SRC_RUNS / rl.source_id(RPM, THETA)
    sp = json.loads((src / "params.json").read_text())
    _, T = pp._t_scales(sp)
    t_rec0 = rl.checkpoint_time(src / "checkpoint.dump") + rl.N_MIX * T
    for tag, (rid, rdir) in CASES.items():
        dump = rg.scratch_dir(rid) / "replay_start.dump"
        assert dump.exists() and (rdir / "flux_rec.0.bin").exists(), (tag, dump, rdir)
        p = rl.replay_params(sp, diag_id(tag), t_rec0, rdir)
        p.update(t_end=1.0 * T, replay_diag=1)
        print(tag, rl._submit(p, dump, "00:30:00"))


def report() -> None:
    for tag in CASES:
        a = np.loadtxt(rg.scratch_dir(diag_id(tag)) / "replay_diag.txt")
        far, allc, um = a[:, 2] / a[:, 4], a[:, 3] / a[:, 4], a[:, 4]
        print(f"{tag:5s} steps {len(a)}  max|uf| {um.max():.4g}  d/umax away from the wall: median {np.median(far):.2e}  max {far.max():.2e}"
              f"   all cells: median {np.median(allc):.2e}  max {allc.max():.2e}")


if __name__ == "__main__":
    {"submit": submit, "report": report}[sys.argv[1]]()
