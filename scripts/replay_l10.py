"""L10 record of ONE flow period (32.5 rpm, theta 7) for the repeated-period test at the finest level (diary 2026-10-22).
Reuses scripts/replay_lf.py (record_params, _submit): restart from the L10 source checkpoint, release after 1 cycle, record one period.
--calibrate: a short job (default 25 min walltime) that only checks that the record binary restores the L10 checkpoint and measures
cost per step and record size per step. Nothing is replayed.
Usage: uv run python scripts/replay_l10.py [--calibrate] [--dry-run] [--walltime HH:MM:SS]
"""
import argparse
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT)]
from scripts import postprocess as pp  # noqa: E402
from scripts import replay_gate as rg  # noqa: E402
from scripts import replay_lf as rl  # noqa: E402

SRC = ROOT / "runs/v_l10_rpm32.5_th7"
REC_ROOT = Path("/oscar/data/dharri15/eaguerov/replay_records/l10")
RID = "lfrec_l10_rpm32.5_th7"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calibrate", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--walltime", default=None)
    ap.add_argument("--replay", action="store_true", help="loop the finished L10 record on the L10 scalar grid (no refinement, no projection)")
    a = ap.parse_args()
    if a.replay:
        sp = json.loads((SRC / "params.json").read_text())
        _, T = pp._t_scales(sp)
        rdir = REC_ROOT / RID
        start = rl.rg.scratch_dir(RID) / "replay_start.dump"
        t_rec0 = rl.checkpoint_time(start)
        lid = "lfrep_l10_rpm32.5_th7"
        p = rl.replay_params(sp, lid, t_rec0, rdir)
        p.update(replay_scalar_level=10, replay_project_flux=0)   # same grid as the flow: the replay is exact (gate F1)
        rl.NTASKS = 32                                            # the record has one flux file per rank (32)
        if a.dry_run:
            print("dry run: would submit", lid, "t_rec0", t_rec0 / T, "cycles", rl.N_REPLAY_CYCLES)
            return
        print("job", rl._submit(p, start, a.walltime or "36:00:00"))
        return
    if not rg.binary_is_current(rl.REC_BIN, rl.SRC_TREE):
        raise SystemExit(f"{rl.REC_BIN} is older than {rl.SRC_TREE}: rebuild first")
    ck = SRC / "checkpoint.dump"
    sp = json.loads((SRC / "params.json").read_text())
    _, T = pp._t_scales(sp)
    t_ck = rl.checkpoint_time(ck)
    print(f"source {SRC.name}: checkpoint at cycle {t_ck / T:.3f}, level params NN={sp.get('NN', sp.get('level'))}")
    rid = RID + ("_cal" if a.calibrate else "")
    rdir = REC_ROOT / rid
    p = rl.record_params(sp, rid, t_ck, rdir)
    rl.NTASKS = 32
    if a.dry_run:
        print("dry run: would submit", rid, "walltime", a.walltime or ("00:25:00" if a.calibrate else "24:00:00"))
        return
    if rdir.exists():
        shutil.rmtree(rdir)
    rdir.mkdir(parents=True)
    print("job", rl._submit(p, ck, a.walltime or ("00:25:00" if a.calibrate else "24:00:00")))


if __name__ == "__main__":
    main()
