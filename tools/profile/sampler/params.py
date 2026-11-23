"""params.json for the sampler check: the profile base case at L6, 2 periods, field snapshots on in period 2.

    params.py OUT_JSON
Diary 2026-11-23 (pre-registration "sampler"): streaming_field.bin (each period end), uv_field.bin (final
checkpoint) and fields/snap_*.bin (period 2) must be bitwise equal between HEAD and the one-reduction sampler.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from tools.profile.params import profile_params  # noqa: E402

if __name__ == "__main__":
    p = profile_params(6, 2) | {"snap_start_cycle": 1.0, "snap_end_cycle": 2.0}
    Path(sys.argv[1]).write_text(json.dumps(p, indent=1))
