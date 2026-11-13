"""Power of the guards: plant one realistic defect in a copy of src/, run the bundle cases it touches, and check
that the guards named for it fail. A guard that passes on its mutant has no power: fix the guard, not the mutant.

    mutants.py run NAME WORK      (on a compute node: tools/bundle/mutants.sbatch NAME WORK)

Each mutant reproduces a defect this project has had or could plausibly get (diary, docs_site/developer/
solver-pitfalls.md). Patches are exact strings: a patch that no longer matches the source fails loudly
(tests/test_bundle_mutants.py), so a mutant cannot silently test nothing.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from tools.bundle import bundle as bb, cases as bc  # noqa: E402


@dataclass(frozen=True)
class Mutant:
    patches: list               # (file in src/, exact old text, new text); old must occur once
    cases: list                 # bundle cases to run (their parents are added)
    must_fail: list             # guard test ids (substring of the pytest node id) that must fail
    flags: dict = field(default_factory=dict)   # extra compile flags per build variant
    note: str = ""


_SOLID = """  bag_solid ();
#if TREE && !REPLAY_FLUX      // a replay keeps the coarse levels it set up in restore_checkpoint()
  restriction ({cs, fs});
#endif
  assert_fs_sane ("post-solid");   // a fully open face must not have fs = 0"""

MUTANTS = {
    # 2026-09-14: face fields are not in the dump; without the recompute ~50 % of open faces are 0 (pitfall fs-after-restore)
    "no_fs_recompute": Mutant([("init.h", _SOLID, "  // (mutant: no cs/fs recompute, no check)")],
                              ["cont", "cold"], ["test_s9_fs_is_sane_after_every_restore[cont]"]),
    # a tracer missing from the dump list: the continuation restarts with no tracer
    "drop_tracer_dump": Mutant([("checkpoint.h", '  write_restart_diagnostic ("pre_dump");\n  dump (file = "checkpoint.dump");',
                                 '  write_restart_diagnostic ("pre_dump");\n  c2.nodump = true;\n  dump (file = "checkpoint.dump");')],
                               ["cont", "cold"], ["test_s10_restore_gives_back_the_dumped_fields", "test_s10_continuation_keeps_the_tracer"]),
    # the replay reads the flux one step late (and the solver's own step check is off, as it once was)
    "replay_off_by_one": Mutant([
        ("replay.h", "    if (rec_pass == 0 && ((nsub == 1 && hi != iter)", "    if (0 && rec_pass == 0 && ((nsub == 1 && hi != iter)"),
        ("replay.h", "    double * buf = (double *) malloc (hn*REC_NPC*sizeof(double));",
         "    static int mutant_skipped = 0;\n"
         "    if (!mutant_skipped) {\n"
         "      mutant_skipped = 1;\n"
         "      fseek (fp_rec, hn*REC_NPC*sizeof(double), SEEK_CUR);\n"
         "      rec_read_head (&hi, &ht, &rec_hdt, &hn);\n"
         "    }\n"
         "    double * buf = (double *) malloc (hn*REC_NPC*sizeof(double));")],
        ["replay", "rec"], ["test_s1_replay_reproduces_the_coupled_scalars"]),
    # the early dump forgets its phase file: the resume takes the nearest boundary (pitfall restart-phase)
    "no_phase_sidecar": Mutant([("checkpoint.h", "  p.nodump = pf.nodump = true;\n  write_phase_sidecar (t);\n  return 1;\n}\n\n// Final checkpoint",
                                 "  p.nodump = pf.nodump = true;\n  return 1;\n}\n\n// Final checkpoint")],
                               ["kill_resume", "cold"], ["test_s5_early_dump_writes_dump_and_phase",
                                                         "test_s5_continuation_reproduces_the_unbroken_run[kill_resume]",
                                                         "test_s12_run_stops_on_the_period_boundary[kill_resume"]),
    # 2026-09-14: refine before solid() on a cross-level restart (pitfall fs-after-restore)
    "prerefine_old": Mutant([], ["xl6", "cold"], ["test_s3_cross_level_warm_start_matches_cold"],
                            flags={"crosslevel": ["-DCROSS_LEVEL_PREREFINE_SOLID=0"]}),
    # cross-case restart without the velocity rescale to the new units
    "no_rescale": Mutant([("init.h", "#if !REPLAY_FLUX\n    rescale_to_new_units ();\n#endif",
                           "#if !REPLAY_FLUX\n    // (mutant: no rescale)\n#endif")],
                         ["caseB_warm", "caseB_cold"], ["test_s4_cross_case_warm_start_matches_cold"],
                         note="S4 is UNSETTLED on the unmutated bundle (diary 2026-11-13): no clean baseline"),
}


def patched_copy(src: Path, dst: Path, patches) -> Path:
    """Copy src to dst and apply exact-string patches; ValueError if an old text does not occur exactly once."""
    shutil.rmtree(dst, ignore_errors=True)
    shutil.copytree(src, dst, ignore=shutil.ignore_patterns("*.o", "BioReactor", "*_preproc*", "_*.c"))
    for f, old, new in patches:
        p = dst / f
        s = p.read_text()
        if s.count(old) != 1:
            raise ValueError(f"{f}: patch text found {s.count(old)} times (must be 1): {old[:60]!r}")
        p.write_text(s.replace(old, new))
    return dst


def verdict(report_lines, must_fail) -> dict:
    """{guard: "failed" | "passed" | "not run"} from pytest -rA summary lines."""
    out = {}
    for g in must_fail:
        hit = [ln for ln in report_lines if g in ln and ln.split(" ")[0] in ("FAILED", "ERROR", "PASSED")]
        if not hit:
            out[g] = "not run"
        else:
            out[g] = "failed" if all(ln.startswith(("FAILED", "ERROR")) for ln in hit) else "passed"
    return out


def caught(v: dict) -> bool:
    return all(x == "failed" for x in v.values())


def run(name: str, work: Path) -> dict:
    m = MUTANTS[name]
    work = work / name
    src = patched_copy(ROOT / "src", work / "src", m.patches)
    builds = {b: f + m.flags.get(b, []) for b, f in bc.BUILDS.items()}
    old = dict(bc.BUILDS)
    bc.BUILDS.clear()
    bc.BUILDS.update(builds)
    try:
        bb.run_bundle(src, work / "bundle", only=m.cases)
    finally:
        bc.BUILDS.clear()
        bc.BUILDS.update(old)
    env = {**os.environ, "MFB_BUNDLE": str(work / "bundle" / "bundle.json"), "MFB_BUNDLE_SRC": str(src)}
    r = subprocess.run([sys.executable, "-m", "pytest", "tests/bundle", "-rA", "-q", "-p", "no:cacheprovider"],
                       cwd=ROOT, env=env, capture_output=True, text=True)
    (work / "pytest.txt").write_text(r.stdout + r.stderr)
    v = verdict(r.stdout.splitlines(), m.must_fail)
    res = {"mutant": name, "guards": v, "caught": caught(v), "note": m.note}
    (work / "verdict.json").write_text(json.dumps(res, indent=1))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run"])
    ap.add_argument("name", choices=sorted(MUTANTS))
    ap.add_argument("work", type=Path)
    a = ap.parse_args()
    res = run(a.name, a.work.resolve())
    print(json.dumps(res, indent=1))
    sys.exit(0 if res["caught"] else 1)


if __name__ == "__main__":
    main()
