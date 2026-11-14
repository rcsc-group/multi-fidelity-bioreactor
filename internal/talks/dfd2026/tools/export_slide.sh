#!/bin/bash
# Export slide(s) of the deck to shots/ as PNG. Usage: tools/export_slide.sh 5   (or 1,4)
# Node 22 is needed by Slidev 53 (the system Node 20 is too old).
export PATH=/oscar/rt/9.6/25/spack/x86_64_v3/node-js-22.16.0-5pf2bfgql2kw7mg3xvibiyip2bwpyf32/bin:$PATH
cd "$(dirname "$0")/.." || exit 1
rm -f shots/*.png
timeout 400 pnpm exec slidev export --format png --range "$1" --output shots 2>&1 | tail -1
ls shots
