"""Visual chaining-fidelity check (diary.md 2026-09-07): a single, 2-segment
chain.py sweep -- segment 0 cold-starts at (32.5rpm, theta=7deg), segment 1
checkpoint-restarts from segment 0's end state into (32.5rpm, theta=4deg),
the same theta-only edge that showed the most dramatic cold-start-beats-
warm-start anomaly in the settling study (L8: cold=2 cycles, warm=58
cycles). Produces a SINGLE continuous set of frames spanning both segments
(chain.py's video path renders each segment separately; frames are
concatenated across segments by scripts/concat_chain_videos.py once both
finish) so the checkpoint transition is visible as a single continuous
video, not two disjoint clips.

Requires BioReactor.c's 2026-09-07 frame-format change (w_bio and
_theta_env_deg appended to each frame) -- rebuild with `make build-video`
before running this.

Usage:
    uv run python scripts/submit_chain_video_demo.py
"""
import sys
from pathlib import Path

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain

chain.validate_params = lambda params: None

cfg = {
    "motion": {"omega_b": 32.5 * 2 * 3.141592653589793 / 60.0, "theta_max": [7.0, 0.0, 0.0]},
    "fidelity": 6,
    "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0},
    "fill_level": 0.5,
    # [FIX 2026-09-07] movies_output() only records frames once t >= t_mix,
    # and t_mix = T_per_st * n_mix_cycles -- the SAME n_mix_cycles used to
    # compute t_end (segment length). With n_mix_cycles=12 that put t_mix
    # almost exactly at t_end, leaving essentially no recording window
    # (segment 0 got 0 frames, segment 1 got 2, diary.md 2026-09-07).
    # Setting n_mix_cycles=0 makes t_mix=0 (records from the start); the
    # full segment length now comes entirely from t_buffer.
    "n_mix_cycles": 0,
    "n_transition_cycles": 0,
    "t_buffer": 9.11,  # ~15 cycles at theta=7's T_per_nd -- same nondim
                        # duration applied to both segments (t_buffer is
                        # shared across the whole chain, chain.py has no
                        # per-segment override)
    "sweep": {"parameter": "theta_max_0", "values": [7.0, 4.0]},
    "mpi": False,
    "videos": True,
    "walltime": "02:00:00",
    "submit": True,
}

job_run_ids = chain.submit_chain(cfg)
print("\nSegments:", job_run_ids)
print(f"run_dirs: runs/{job_run_ids[0][0]}/  (cold start)")
print(f"          runs/{job_run_ids[1][0]}/  (warm start, checkpointed from above)")
