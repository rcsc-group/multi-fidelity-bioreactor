"""Rows per level for kLa_1T_10 and tau_95 (after FAILED-GATE removal), in the campaign domain."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import campaign as C  # noqa: E402
g = C.gate_all_runs()
for q in ("kLa_1T_10", "tau_95"):
    d = C.exclude_failed_gate(C.load_db(q), g)
    print(q, d.groupby("level").size().to_dict())
