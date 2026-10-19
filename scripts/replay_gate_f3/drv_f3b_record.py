"""F3b step 1: L8 coupled run (T5 protocol, release at cycle 80) that records 22 post-release cycles at 16 ranks.
Record in /oscar/data/dharri15/eaguerov/replay_records/f3b_l8_rpm32.5/ ; node2356 excluded."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts import replay_gate as rg

REC = Path("/oscar/data/dharri15/eaguerov/replay_records/f3b_l8_rpm32.5")
REC.mkdir(parents=True, exist_ok=True)
p = rg.make_params("f3b_coupled", rg.BIN_DIR / "BioReactor-recL", fidelity=8, n_mix_cycles=80, post_cycles=22,
                   n_rec_periods=22, rec_dir=str(REC))
# T5 geometry/forcing are the make_params defaults (a=0.25, b=0.03575, n=8, theta 7, 32.5 rpm)
job = rg.submit(p, walltime="08:00:00", ntasks=16, mem="4G", exclude="node2356")
print("job", job, rg.wait(job), flush=True)
