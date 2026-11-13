"""Repo audit for the lean-repo spec: per tracked file, last commit date and whether any OTHER tracked file names it.
Output: TSV (path, last_date, n_referrers, lines) + a summary by top-level folder. Read-only.
Usage: python repo_audit.py REPO OUT.tsv
"""
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

repo, out = Path(sys.argv[1]), Path(sys.argv[2])
git = lambda *a: subprocess.run(["git", "-C", str(repo), *a], capture_output=True, text=True).stdout  # noqa: E731

files = git("ls-files").split("\n")
files = [f for f in files if f]

# last commit date per file: one pass over the whole log (newest first; first sighting = last change)
last = {}
date = None
for line in git("log", "--name-only", "--format=@@%cs").split("\n"):
    if line.startswith("@@"):
        date = line[2:]
    elif line and line not in last:
        last[line] = date

# referrers: does any other tracked text file mention the file's stem (module import) or its name?
TEXT = (".py", ".sh", ".sbatch", ".md", ".toml", ".yml", ".yaml", ".c", ".h", ".json", ".cfg", ".txt")
corpus = {}
for f in files:
    if f.endswith(TEXT) and not f.startswith("experiments/"):
        try:
            corpus[f] = (repo / f).read_text(errors="ignore")
        except OSError:
            pass
rows = []
for f in files:
    p = Path(f)
    keys = {p.name}
    if p.suffix == ".py":
        keys.add(p.stem)
    n = sum(1 for g, txt in corpus.items() if g != f and any(k in txt for k in keys))
    try:
        nl = sum(1 for _ in open(repo / f, errors="ignore")) if f.endswith(TEXT) else 0
    except OSError:
        nl = 0
    rows.append((f, last.get(f, "?"), n, nl))
with open(out, "w") as fh:
    fh.write("path\tlast\tn_ref\tlines\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")

agg = defaultdict(lambda: [0, 0, 0, 0])
for f, d, n, nl in rows:
    top = f.split("/")[0] if "/" in f else "(root)"
    a = agg[top]; a[0] += 1; a[1] += nl; a[2] += (n == 0); a[3] += (d < "2026-10-13")
print("folder files lines unreferenced older_than_30d")
for k, (c, nl, u, o) in sorted(agg.items(), key=lambda kv: -kv[1][0]):
    print(f"{k:14s} {c:5d} {nl:7d} {u:5d} {o:5d}")
