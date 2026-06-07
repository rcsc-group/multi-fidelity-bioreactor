"""Render and commit the canonical Kim et al. 2024 bioreactor case.

Waits for canonical_kim2024 run to complete, renders the VOF video,
copies all outputs to docs/canonical_case/, and prints the git commands
to commit.

Usage:
    uv run python scripts/commit_canonical.py
"""
import json, shutil, subprocess, time
from pathlib import Path

PROJECT_ROOT = Path(__file__).parents[1]
RUN_DIR      = PROJECT_ROOT / "runs" / "canonical_kim2024"
SCRATCH_DIR  = Path("/oscar/scratch/eaguerov/mpi_runs/canonical_kim2024")
DOCS_DIR     = PROJECT_ROOT / "docs" / "canonical_case"
DOCS_DIR.mkdir(parents=True, exist_ok=True)

print("Waiting for canonical_kim2024 to complete...")
deadline = time.monotonic() + 7200
while time.monotonic() < deadline:
    if (RUN_DIR / "results.json").exists():
        break
    time.sleep(30)
    print(f"  still running... ({int((deadline - time.monotonic())/60)} min left)")
else:
    raise TimeoutError("canonical_kim2024 did not complete within 2 hours")

print("Run complete. Rendering video...")

# Render from scratch dir (frames are there before rsync)
src = SCRATCH_DIR if SCRATCH_DIR.exists() else RUN_DIR
subprocess.run(
    ["uv", "run", "python", "scripts/render_videos.py", str(src)],
    cwd=PROJECT_ROOT, check=True
)

# Copy outputs to docs/canonical_case/
for fname in ["volume_fraction.mp4", "volume_fraction_lab.mp4",
              "results.json", "params.json",
              "logstats.dat", "tr_oxy.dat", "normf.dat"]:
    for src_dir in [RUN_DIR, SCRATCH_DIR]:
        src_f = src_dir / fname
        if src_f.exists():
            shutil.copy2(src_f, DOCS_DIR / fname)
            print(f"  copied {fname}")
            break

# Print KPIs
res = json.loads((DOCS_DIR / "results.json").read_text())
print("\nCanonical case KPIs (omega_b=3.93, theta=7°, fill=0.5, fidelity=7):")
for k in ["kLa_25", "kLa_50", "vel_rms_qss", "kla_fit_rmse_25", "dtmix_0.50"]:
    print(f"  {k}: {res.get(k)}")

print("\nTo commit, run:")
print("  git add docs/canonical_case/")
print("  git commit -m 'docs: canonical Kim et al. 2024 case video and KPIs'")
