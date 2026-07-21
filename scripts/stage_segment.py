"""Safely stage a chain-restart segment's checkpoint + params into scratch.

Root cause this guards against: checkpoint.dump existing on disk is NOT proof
a segment finished. A job killed seconds after starting a restart still has
a checkpoint.dump on disk -- it's just an untouched copy of the segment's own
INPUT, not its output. The only reliable proof a segment ran to completion is
a well-formed results.json, written exclusively by postprocess.py at the end
of a real run. Every restart must go through this module rather than an
ad-hoc `cp`, so that guarantee is never bypassed by hand.

Usage:
    uv run python scripts/stage_segment.py <prev_run_id> <next_run_id> \
        [--runs-root runs] [--scratch-root /oscar/scratch/eaguerov/mpi_runs]
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

DEFAULT_SCRATCH_ROOT = Path("/oscar/scratch/eaguerov/mpi_runs")


class SegmentNotCompleteError(RuntimeError):
    """Raised when a segment claimed as a restart source has no valid results.json."""


def _verify_segment_complete(run_dir: Path) -> None:
    results_file = run_dir / "results.json"
    if not results_file.exists():
        raise SegmentNotCompleteError(
            f"{run_dir}: no results.json -- this segment has not actually "
            f"finished, even if checkpoint.dump exists (a checkpoint.dump can "
            f"be present from an unrun/killed restart job; only results.json, "
            f"written by postprocess.py, proves completion)"
        )
    try:
        results = json.loads(results_file.read_text())
    except json.JSONDecodeError as e:
        raise SegmentNotCompleteError(f"{run_dir}: results.json is malformed: {e}") from e
    if not isinstance(results, dict) or not results:
        raise SegmentNotCompleteError(f"{run_dir}: results.json is empty/invalid")


def stage_segment(
    prev_run_id: str,
    next_run_id: str,
    *,
    runs_root: Path,
    scratch_root: Path = DEFAULT_SCRATCH_ROOT,
) -> Path:
    """Stage next_run_id's params + prev_run_id's checkpoint into scratch.

    Returns the scratch directory staged into. Raises SegmentNotCompleteError
    if prev_run_id has no valid results.json, or FileNotFoundError if the
    required source files don't exist.
    """
    prev_dir = runs_root / prev_run_id
    next_dir = runs_root / next_run_id

    _verify_segment_complete(prev_dir)

    prev_checkpoint = prev_dir / "checkpoint.dump"
    if not prev_checkpoint.exists():
        raise FileNotFoundError(f"{prev_checkpoint}: prev segment's checkpoint.dump not found")

    next_params = next_dir / "params.json"
    if not next_params.exists():
        raise FileNotFoundError(f"{next_params}: next segment's params.json not found")

    staged_dir = scratch_root / next_run_id
    staged_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(prev_checkpoint, staged_dir / "checkpoint.dump")
    shutil.copy2(next_params, staged_dir / "params.json")
    return staged_dir


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prev_run_id")
    parser.add_argument("next_run_id")
    parser.add_argument("--runs-root", type=Path, default=Path(__file__).parents[1] / "runs")
    parser.add_argument("--scratch-root", type=Path, default=DEFAULT_SCRATCH_ROOT)
    args = parser.parse_args()

    try:
        staged_dir = stage_segment(
            args.prev_run_id, args.next_run_id,
            runs_root=args.runs_root, scratch_root=args.scratch_root,
        )
    except (SegmentNotCompleteError, FileNotFoundError) as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Staged {args.next_run_id} <- {args.prev_run_id} at {staged_dir}")


if __name__ == "__main__":
    main()
