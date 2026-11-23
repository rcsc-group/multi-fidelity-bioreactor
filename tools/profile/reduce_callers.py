"""Count MPI reductions by caller: patch qcc's generated _BioReactor.c, then map the counts back to source lines.

    reduce_callers.py patch _BioReactor.c          (in place; before mpicc)
    reduce_callers.py report BINARY COUNTS STEPS   (after the run; addr2line on the caller addresses)

The patch makes mpi_all_reduce_array (every foreach reduction goes through it) record its two return addresses.
Rank 0 rewrites reduce_callers.txt every 20000 calls and at exit, so a run killed by its time limit still
leaves counts. Basilisk itself is not touched: only the generated copy is. Diary 2026-11-23 (profiling).
"""
import subprocess
import sys
from collections import Counter
from pathlib import Path

HOOK = r"""
#define RC_N 1024
static void * rc_a0[RC_N], * rc_a1[RC_N];
static long rc_cnt[RC_N], rc_total = 0;
static int rc_used = 0, rc_rank = -1;
static void rc_dump (void) {
  if (rc_rank != 0) return;
  FILE * f = fopen ("reduce_callers.txt", "w");
  if (!f) return;
  for (int i = 0; i < rc_used; i++)
    fprintf (f, "%p %p %ld\n", rc_a0[i], rc_a1[i], rc_cnt[i]);
  fclose (f);
}
static void rc_hit (void * a0, void * a1) {
  if (rc_rank < 0) { MPI_Comm_rank (MPI_COMM_WORLD, &rc_rank); atexit (rc_dump); }
  int i;
  for (i = 0; i < rc_used; i++)
    if (rc_a0[i] == a0 && rc_a1[i] == a1) break;
  if (i == rc_used && rc_used < RC_N) { rc_a0[i] = a0; rc_a1[i] = a1; rc_cnt[i] = 0; rc_used++; }
  if (i < rc_used) rc_cnt[i]++;
  if (++rc_total % 20000 == 0) rc_dump ();
}
"""
SIG = "void mpi_all_reduce_array (void * v, MPI_Datatype datatype, MPI_Op op, int elem)\n{"


def patch(path: Path) -> None:
    s = path.read_text()
    if s.count(SIG) != 1:
        raise ValueError(f"signature found {s.count(SIG)} times (must be 1)")
    s = s.replace(SIG, HOOK + SIG + "\n  rc_hit (__builtin_return_address (0), __builtin_return_address (1));")
    path.write_text(s)


def where(binary: str, addrs: list) -> dict:
    """{addr: 'func file:line'} via addr2line (non-PIE binary, -g)."""
    out = subprocess.run(["addr2line", "-f", "-C", "-e", binary] + addrs, capture_output=True, text=True).stdout
    ln = out.splitlines()
    return {a: f"{ln[2 * k]} {Path(ln[2 * k + 1].split(' ')[0]).name}:{ln[2 * k + 1].split(':')[-1].split(' ')[0]}"
            for k, a in enumerate(addrs)}


def report(binary: str, counts: Path, steps: int) -> None:
    rows = [ln.split() for ln in counts.read_text().splitlines() if ln.strip()]
    loc = where(binary, sorted({r[0] for r in rows} | {r[1] for r in rows}))
    by = Counter()
    for a0, a1, n in rows:
        by[(loc[a0], loc[a1])] += int(n)
    total = sum(by.values())
    print(f"total {total} reductions on rank 0, {steps} steps -> {total / steps:.1f} per step")
    for (l0, l1), n in by.most_common(40):
        print(f"{n / steps:8.2f}/step  {100 * n / total:5.1f} %  {l0}  <-  {l1}")


if __name__ == "__main__":
    if sys.argv[1] == "patch":
        patch(Path(sys.argv[2]))
    elif sys.argv[1] == "report":
        report(sys.argv[2], Path(sys.argv[3]), int(sys.argv[4]))
    else:
        sys.exit("usage: reduce_callers.py patch FILE | report BINARY COUNTS STEPS")
