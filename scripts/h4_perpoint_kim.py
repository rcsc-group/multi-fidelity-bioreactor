"""Per-point error of the chain at the theta 5.5 test set, Kim statistics (reads test_h4_testset_kim.json). Usage: uv run python scripts/h4_perpoint_kim.py"""
import json
import sys

sys.path[:0] = ["."]
from scripts.hydro_dataset_v2 import qois  # noqa: E402

d = json.load(open("experiments/multifidelity/test_h4_testset_kim.json"))
rpms = [17.5, 20, 22.5, 27.5, 30, 32.5]
obs = {r: qois(f"x8_l10_rpm{r:g}_th5.5") for r in rpms}
for q in ("tau_amp_t", "edr_max_t"):
    print(q, " ".join(f"{r:g}rpm {100 * (p / obs[r][q] - 1):+.0f}%" for r, p in zip(rpms, d[q]["A"]["pred"])))
