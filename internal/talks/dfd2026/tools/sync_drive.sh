#!/bin/bash
# Push the deck source and its rendered PDF to Drive: "Elvis/Presentations/APS DFD/2026/".
# Default is a DRY RUN that lists what would change. Run with --go only after Elvis approves.
# rclone copy never deletes on Drive; files there that are not here are kept.
set -euo pipefail
module load rclone/1.70.2-mfg5
cd "$(dirname "$0")/.."
DST="gdrive:Elvis/Presentations/APS DFD/2026/"
FLAGS=(--include "slides.md" --include "style.css" --include "*.pdf" --max-depth 1)
if [[ "${1:-}" == "--go" ]]; then
  rclone copy . "$DST" "${FLAGS[@]}" -v
else
  rclone copy . "$DST" "${FLAGS[@]}" --dry-run 2>&1 | grep -E "NOTICE|Skipped" || true
  echo "dry run only; rerun with --go to copy"
fi
