"""Does the mixing-time bias change with the condition? (pre-registered in diary.md, 2026-10-29)
Same L8 DNS per condition, two binaries: OLD = production before the fix (seed f > 0.5, dtmix from c2*f),
NEW = fixed production (Kim's seed f == 1, dtmix from bulk-liquid sums with sigma^2_max at the release).
  build             BioReactor-mixold (-DSEED_KIM=0) and BioReactor-mixnew (defaults), MPI
  smoke             L6 NEW run at 20 rpm theta 7, 3 cycles after release (pipeline check, not the experiment)
  submit            the six L8 runs (3 conditions x OLD/NEW), 200 s of physical time after the release
  report            dtmix of every run and the ratio NEW/OLD per condition; writes docs/meetings/prereg/mixfix.txt
Usage: uv run python scripts/mixfix_check.py build | smoke | submit | report
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import postprocess as pp  # noqa: E402
from scripts import prereg_offline as po  # noqa: E402
from scripts import replay_gate as rg  # noqa: E402

BIN = {"old": rg.BIN_DIR / "BioReactor-mixold", "new": rg.BIN_DIR / "BioReactor-mixnew"}
FLAGS = {"old": ["-DSEED_KIM=0"], "new": []}
CONDS = [(17.5, 5.5), (25.0, 7.0), (32.5, 7.0)]      # (rpm, theta): span of the sweep / test set
N_MIX = 20
T_AFTER_S = 200.0                                    # physical seconds after the release
NT = 16
OUT = ROOT / "docs/meetings/prereg"


def rid(kind: str, rpm: float, theta: float, level: int = 8) -> str:
    return f"mix{level}_{kind}_{rpm:g}_{theta:g}".replace(".", "p")


def params(kind: str, rpm: float, theta: float, level: int, post_cycles: float | None = None) -> dict:
    p = rg.make_params(rid(kind, rpm, theta, level), BIN[kind], fidelity=level, n_mix_cycles=N_MIX,
                       post_cycles=1.0, rpm=rpm, theta=theta)
    T_bio, t_per = pp._t_scales(p)
    cyc = post_cycles if post_cycles is not None else math.ceil(T_AFTER_S / (t_per * T_bio))
    p["t_end"] = (N_MIX + cyc) * t_per
    return p


def build() -> None:
    for k, b in BIN.items():
        if b.exists():
            raise RuntimeError(f"{b} exists: never overwrite a binary; remove it by hand if it is yours")
        print(rg.build_binary(b.name.removeprefix("BioReactor-"), FLAGS[k], mpi=True))


def smoke() -> None:
    print(po._submit(params("new", 20.0, 7.0, 6, post_cycles=3.0), NT, "01:00:00"))


def submit() -> None:
    jobs = {}
    for rpm, theta in CONDS:
        for kind in ("old", "new"):
            p = params(kind, rpm, theta, 8)
            jobs[p["run_id"]] = po._submit(p, NT, "12:00:00")
    (rg.WORK / "mixfix_jobs.json").write_text(json.dumps(jobs))
    print(jobs)


def report() -> None:
    rows = []
    for rpm, theta in CONDS:
        d = {}
        for kind in ("old", "new"):
            run = rg.RUNS / rid(kind, rpm, theta)
            sp = json.loads((run / "params.json").read_text())
            d[kind] = pp._compute_mixing_metrics(run, sp, use_bulk=(kind == "new"))   # old = old method
        r = {k: d["new"][k] / d["old"][k] for k in ("dtmix_0.50", "dtmix_0.75", "dtmix_0.95")}
        rows.append(f"{rpm:5g} rpm theta {theta:g}: "
                    + "  ".join(f"{k[-4:]} old {d['old'][k]:7.2f} s new {d['new'][k]:7.2f} s ratio {r[k]:.3f}"
                                for k in r))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "mixfix.txt").write_text("\n".join(rows) + "\n")
    print("\n".join(rows))


if __name__ == "__main__":
    {"build": build, "smoke": smoke, "submit": submit, "report": report}[sys.argv[1]]()
