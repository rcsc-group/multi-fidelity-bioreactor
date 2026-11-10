"""Paired L7 -> L8 runs at the 12 good off-line L9 sites (bo1, N_WARM = 30) (user go, 2026-11-10 ~16:40).

Why: off theta = 7 the L8 -> L9 link for kLa is UNMEASURED. The only off-line sites with both L8 and L9 are the c2d runs,
dropped (released after 12 settle cycles; the L9 wave needs ~27). The bo1 L9 runs have no L8 at the same site.
Protocol = the live chain arm (bo4): L7 cold, then L8 warm from that L7, same forcing. Run ids pr8_l7_*, pr8_l8_*.

PRE-REGISTERED (written before any result; revised 16:50 before any run started, user: "Pearson is a disservice to
mfbml if the transfer is not linear"): the score is the one mfbml itself uses - prediction of held-out L9.
Leave-one-site-out over the 12 sites (+ the bo4 chain sites with L8 and L9, reported separately), L8 value known at
the held-out site (nested design): RMSE (ln) of (a) L9-only GP, (b) mfbml T linear in f8, (c) T quadratic in f8,
(d) T linear in f8 + theta; all with the GP discrepancy Delta(x). L8 "carries L9 information off-line" if the best
of (b)-(d) has RMSE <= 0.8 x (a) for kLa_1T_10. Form chosen by inner LOO only (no peeking at the outer score).
Secondary, descriptive only: Spearman rho (monotone, not linear). Same scores for tau_95 and dtmix_0.25.
Power: n = 12, so a 20 % RMSE gap is near the draw scatter; this is a screen, not a proof.

Driver: own task file state/tasks_pair.json, reuses loop.advance (the running autopilot's tasks.json is not touched).
  uv run python scripts/bo2d/pair_l8.py init|advance [--go]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.bo2d import loop as L  # noqa: E402

SITES = [(15.45, 5.0), (17.025, 4.8), (17.925, 3.15), (18.6, 4.6), (20.175, 4.4), (21.525, 3.15),
         (30.75, 5.5), (33.225, 4.5), (33.45, 3.05), (34.35, 2.7), (35.475, 2.35), (36.375, 2.0)]
TAG = "pr8"
L.TASKS = L.S.STATE / "tasks_pair.json"
L.TAGS["pair"] = TAG
_src = L._source_for


def _source_for(t, step):
    if t["arm"] == "pair" and step == "l8":
        return L.C.run_id(7, *t["x"], tag=TAG)
    return _src(t, step)


L._source_for = _source_for

if __name__ == "__main__":
    go = "--go" in sys.argv
    if sys.argv[1] == "init":
        assert not L.TASKS.exists(), "tasks_pair.json exists"
        L.save_tasks([dict(arm="pair", obj="pair", x=[r, t], steps=["l7", "l8"], i=0, jobs=[], runs=[], state="new")
                      for r, t in SITES])
        print("wrote", L.TASKS, len(SITES), "tasks")
    else:
        L.advance(go)
        left = [t for t in L.load_tasks() if t["state"] not in ("done", "failed")]
        print("left", len(left))
