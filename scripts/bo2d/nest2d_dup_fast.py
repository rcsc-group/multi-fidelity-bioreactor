"""Speed test by duplication (2026-11-11 ~19:00, user: "how do you know it is the node's fault? ... send in parallel on
a fast node, if it is truly faster, then cancel the slow one").
H: ns_l9_rpm15.42_th4.66 (job 8868642, node1818 = Intel Cascade Lake, 48 cores) runs at ~880 s per time unit because of
the node; the other nested L9 ran at 250-420 s per unit on AMD Genoa nodes. Test: an identical warm start from the same
L8 parent, tag nsf, restricted to Genoa nodes (scontrol Features=genoa after submit). Falsified if its s/unit is not
clearly below ~880 after the same physical-time span. Keep the faster one; the loser is cancelled and its dir renamed.
  uv run python scripts/bo2d/nest2d_dup_fast.py --go
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.bo2d import submit as S  # noqa: E402

print(S.warm(9, "ns_l8_rpm15.42_th4.66", 15.423, 4.655, arm="nest", go="--go" in sys.argv, tag="nsf"))
