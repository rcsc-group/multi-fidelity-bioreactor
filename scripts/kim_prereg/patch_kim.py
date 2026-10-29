"""Patch the scratch copy of Kim's code (kim_pipeline/prereg) for the offline-vs-DNS test (diary 2026-10-29).
Every change is behind a switch that defaults to 0 = his code unchanged. Writes prereg_src/changes.diff.
Usage: uv run python /oscar/scratch/eaguerov/kim_pipeline/prereg_src/patch_kim.py
"""
from __future__ import annotations

import difflib
import shutil
from pathlib import Path

ROOT = Path("/oscar/scratch/eaguerov/kim_pipeline/prereg")
REF = Path("/oscar/scratch/eaguerov/kim_pipeline/repo_pristine_ref")
HDR = '#include "/oscar/scratch/eaguerov/kim_pipeline/prereg_src/prereg_MK.h"'


def once(text: str, old: str, new: str) -> str:
    n = text.count(old)
    assert n == 1, f"anchor occurs {n} times: {old[:80]!r}"
    return text.replace(old, new)


# producer: two switches the case file can set (BR_PREREG_REC, BR_PREREG_DIAG), header at the end
p = ROOT / "codes_bio_run/BioReactor.c"
s = (REF / "codes_bio_run/BioReactor.c").read_text()
s = once(s, "#define OUT_DUMP_BINARY   1\n",
         "#define OUT_DUMP_BINARY   1\n#define PREREG_REC        0   // offline-vs-DNS test: record u, f, uf each step\n"
         "#define PREREG_DIAG       0   // offline-vs-DNS test: per-step scalar sums\n")
s = s.rstrip("\n") + "\n\n" + HDR + "\n"
p.write_text(s)

# offline: mark the file, disable restore_all + Field_define and logstats in the test modes, header at the end
q = ROOT / "codes/TracerTest_test.c"
t = (REF / "codes/TracerTest_test.c").read_text()
t = once(t, "#ifndef NFILES\n", "#define PREREG_OFFLINE 1\n#ifndef PREREG_MODE\n#define PREREG_MODE 0\n#endif\n"
         "#ifndef PREREG_DIAG\n#define PREREG_DIAG 0\n#endif\n#ifndef NFILES\n")
t = once(t, "// event restore_all(i++){\nevent restore_all(",
         "#if !PREREG_MODE\n// event restore_all(i++){\nevent restore_all(")
t = once(t, "  dt = dtnext(timestep (uf, DT));\n  TM_STOP (TM_FIELDDEF);\n}\n",
         "  dt = dtnext(timestep (uf, DT));\n  TM_STOP (TM_FIELDDEF);\n}\n#endif // !PREREG_MODE\n")
t = once(t, "//  Log performance and runtime\nevent logstats (t+=0.1; t <= t_end_file) {",
         "#if !PREREG_MODE\n//  Log performance and runtime\nevent logstats (t+=0.1; t <= t_end_file) {")
t = once(t, "// Output files\n#if OUT_FILES\nevent out_files(", "#endif // !PREREG_MODE\n\n// Output files\n#if OUT_FILES\nevent out_files(")
t = once(t, "#if RESTORE\nevent initial_conditions (i = 0) {",
         "#if RESTORE\n#if !PREREG_MODE   // test modes seed in prereg_MK.h, where the producer seeds\n"
         "event initial_conditions (i = 0) {")
t = once(t, "  boundary({c, c1, c2, c3, oxy});\n}\n\n#if OXYGEN_REPLE",
         "  boundary({c, c1, c2, c3, oxy});\n}\n#endif // !PREREG_MODE\n\n#if OXYGEN_REPLE")
# qcc does not resolve the absolute include from the offline build dirs: local include + a copy in codes/
t = t.rstrip("\n") + "\n\n" + '#include "prereg_MK.h"' + "\n"
shutil.copy2(Path(__file__).parent / "prereg_MK.h", ROOT / "codes/prereg_MK.h")
q.write_text(t)

diff = []
for rel in ("codes_bio_run/BioReactor.c", "codes/TracerTest_test.c"):
    a = (REF / rel).read_text().splitlines(keepends=True)
    b = (ROOT / rel).read_text().splitlines(keepends=True)
    diff += difflib.unified_diff(a, b, f"a/{rel}", f"b/{rel}")
(Path(__file__).parent / "changes.diff").write_text("".join(diff))
print("patched; diff lines:", len(diff))
