"""Task 3 (2026-11-04): re-measure the 2026-10-21 K-test (22 cycles, 32.5 rpm, L8) with the CURRENT functions, read-only.
Usage: uv run python scripts/check_replay_kla25_f3b.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import replay_gate as rg  # noqa: E402
from scripts.check_replay_kla25_metrics import current, QOIS  # noqa: E402

ref_dir = rg.run_dir("f3b_coupled")
ref = current(ref_dir)
print("coupled (current):", {q: round(ref[q], 3) for q in QOIS}, "span", round(ref["span_cycles"], 1))
for k in (1, 2, 4):
    d = rg.run_dir(f"f3b_loop{k}")
    r = current(d)
    dev = {q: 100 * (r[q] / ref[q] - 1) for q in QOIS}
    print(f"K={k} (current):", "  ".join(f"{q} {dev[q]:+.1f}" for q in QOIS), "span", round(r["span_cycles"], 1))
