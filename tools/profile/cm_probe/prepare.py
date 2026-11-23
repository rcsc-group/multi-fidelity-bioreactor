"""Make the source copy for one cm-probe arm: src/ + cm_probe.h appended to BioReactor.c.

    prepare.py ARM DST       ARM = keep (HEAD) | drop (restriction({cm}) outside the gate removed)
                             | once / every (drop + restriction({cm}) after fractions_cleanup at step 0 / each step)
Diary 2026-11-23 (profiling, pre-registration "cm restriction").
"""
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from tools.bundle.mutants import patched_copy  # noqa: E402

DROP = ("henry_oxy2.h", "    restriction ({cm});\n    if (cabs > 0.) {", "    if (cabs > 0.) {")


def prepare(arm: str, dst: Path) -> Path:
    patches = {"keep": [], "drop": [DROP], "once": [DROP], "every": [DROP]}[arm]
    patched_copy(ROOT / "src", dst, patches)
    extra = ["cm_probe.h"] + (["cm_after_cleanup.h"] if arm in ("once", "every") else [])
    main = dst / "BioReactor.c"
    text = main.read_text() + ("\n#define CM_ONCE 1\n" if arm == "once" else "\n")
    for h in extra:
        shutil.copy(Path(__file__).with_name(h), dst / h)
        text += f'#include "{h}"   // profiling probe (tools/profile/cm_probe)\n'
    main.write_text(text)
    return dst


if __name__ == "__main__":
    prepare(sys.argv[1], Path(sys.argv[2]))
