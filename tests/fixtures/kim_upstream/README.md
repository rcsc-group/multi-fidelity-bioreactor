# Vendored Kim et al. upstream driver code

These files are Minki Kim et al.'s own published driver code
(https://github.com/rcsc-group/BioReactor/tree/main/DriverCodes), fetched
verbatim and patched with the **minimal, individually-justified changes**
needed to compile against this project's 2026 Basilisk install. Every
change is marked inline with a `[PROJECT ADDED/CHANGED]` comment; nothing
else was touched. See `diary.md` (repo root) for the full derivation —
search for `kim_upstream_clean`, the scratch directory this fixture was
copied from.

The four changes, all in `BioReactor.c` / `henry_oxy2.h`:

1. `t_end` truncated 250 → 24 (test-only truncation, no physics effect —
   just how long the run continues past the comparison window).
2. `L0 = 1.[0]; DT = HUGE[0];` — dimensional-analysis annotations current
   Basilisk requires on literals; upstream's own values (`LL`, unset
   default) unchanged.
3. `henry_oxy2.h`: `set_prolongation`/`set_restriction` API rename —
   current Basilisk's embed.h dropped the old `.dirty=true` /
   `s.refine=s.prolongation=...` API this file used.
4. `henry_oxy2.h`: `q.embed_flux = NULL;` — a real, independent bug in
   Kim et al.'s own code (uninitialized struct field read later via
   `if (!q.embed_flux && ...)`); not a Basilisk-version issue.

**Why this exists:** `tests/verification/test_kim_upstream_comparison.py`
uses this fixture to report (as a non-fatal warning, never a test
failure) how far this project's own fork has drifted from a
close-to-literal reproduction of Kim et al.'s published code, for a
condition and resolution matched as closely as possible to this fork's
own regression tests. This is NOT a correctness gate — this project's
fork intentionally differs from upstream in several documented ways
(ramp duration, tank cross-section shape, liquid-volume convention; see
`diary.md`, 2026-07-30 entries), so a discrepancy here is expected, not a
bug. It exists purely so a large, unexpected *change* in that discrepancy
(e.g. from a future refactor) is visible to whoever runs it, instead of
going unnoticed. This test is `medium`-marked and runs automatically in
GitHub Actions CI (`medium-tests` job) — Basilisk is built from source
there already, so no OSCAR access is needed.

**Do not edit these files to "fix" the discrepancy against our own fork.**
If Kim's own code needs a change to keep compiling against a newer
Basilisk, make it here as its own documented, minimal, individually
justified diff — same discipline as the rest of this fixture — and record
it in this README and in `diary.md`.
