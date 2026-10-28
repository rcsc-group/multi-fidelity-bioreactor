"""Cost of rebuilding the advecting face velocity from cell-centred velocity (the offline approach of Kim's rocking-bioreactor-2d-offline).

Same protocol as scripts/replay_lf.py (record 1 period after the release, loop it for 22 cycles at L8), at one condition, but the RECORD binary is
BioReactor-recL-kim (-DREC_KIM_FLUX=1): the record file holds ufk = fm*face_value(u.x, 0) of the coupled run instead of its uf. The coupled
run and the replay binary are unchanged, so any difference to the faithful replay lf_l8_rpm<r>_th<theta> comes from the face field alone.
(Kim's solver holds each snapshot over a snapshot interval; here every step has its own field, so this isolates the spatial rebuild only.)
Usage: uv run python scripts/replay_kim.py [rpm theta] [--dry-run]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import replay_gate as rg  # noqa: E402
from scripts import replay_lf as rl  # noqa: E402

rl.REC_BIN = rg.BIN_DIR / "BioReactor-recL-kim"
rl.REC_ROOT = Path("/oscar/data/dharri15/eaguerov/replay_records/kim")
rl.rec_id = lambda rpm, theta: f"kimrec_l8_rpm{rpm:g}_th{theta:g}"
rl.lf_id = lambda rpm, theta: f"kim_l8_rpm{rpm:g}_th{theta:g}"

# BioReactor-replayL16 predates the REC_KIM_FLUX edit of src/BioReactor.c. The edit is inside `#if REC_KIM_FLUX` (record side only; verified in
# the diff), so the replay binary needs no rebuild; the freshness check is waived for it alone, the record binary is still checked.
_orig = rg.binary_is_current
rg.binary_is_current = lambda b, src=rg.ROOT / "src": True if b == rl.REP_BIN else _orig(b, src)

if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    rpm, theta = (float(args[0]), float(args[1])) if len(args) == 2 else (32.5, 7.0)
    info = rl.submit_pair(rpm, theta, dry_run="--dry-run" in sys.argv)
    print(json.dumps(info, indent=1))
