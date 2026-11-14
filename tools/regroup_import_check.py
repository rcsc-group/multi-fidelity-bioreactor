"""Check 3 of the 2026-11 script regroup: every moved script's imports resolve from its new location wherever they
resolved from its old one.

For each moved script, run only its top-level import statements and sys.path set-up (plus the assignments those use)
twice in a fresh interpreter: the old text (git show REV:scripts/NAME.py) as if at scripts/NAME.py, and the new
file at scripts/<dir>/NAME.py. sys.path[0] is the script's own folder, as for `python scripts/.../NAME.py`.
A regression is "imports worked at the old place, fail at the new place".

Usage: python tools/regroup_import_check.py MOVES_JSON REV [N_WORKERS] [--also FILE ...]   (repo root, compute node)
  --also: files that did not move but whose imports were rewritten; checked old text vs new text in place.
"""
import ast
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path.cwd()

RUNNER = r'''
import sys, runpy
src, fake = sys.argv[1], sys.argv[2]
sys.path[0] = str(__import__("pathlib").Path(fake).parent)
sys.argv = [fake]
g = {"__file__": fake, "__name__": "__regroup_check__"}
exec(compile(open(src).read(), fake, "exec"), g)
'''


def _names_used(node) -> set:
    return {n.id for n in ast.walk(node) if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)}


def _names_set(node) -> set:
    return {n.id for t in getattr(node, "targets", [getattr(node, "target", None)]) if t is not None
            for n in ast.walk(t) if isinstance(n, ast.Name)}


def imports_only(text: str) -> str:
    """Top-level imports and sys.path set-up, plus (closure) every top-level assignment they use."""
    tree = ast.parse(text)
    keep = set()
    for i, node in enumerate(tree.body):
        s = ast.unparse(node)
        if isinstance(node, (ast.Import, ast.ImportFrom)) or "sys.path" in s:
            keep.add(i)
        elif isinstance(node, ast.Try) and any(isinstance(b, (ast.Import, ast.ImportFrom)) for b in node.body):
            keep.add(i)
    assigns = {i: n for i, n in enumerate(tree.body) if isinstance(n, (ast.Assign, ast.AnnAssign))}
    while True:
        need = set().union(*(_names_used(tree.body[i]) for i in keep)) if keep else set()
        add = {i for i, n in assigns.items() if i not in keep and _names_set(n) & need}
        if not add:
            break
        keep |= add
    return ast.unparse(ast.Module(body=[tree.body[i] for i in sorted(keep)], type_ignores=[]))


def run(text: str, fake: Path, tag: str, name: str, tmp: Path):
    src = tmp / f"{name}.{tag}.py"
    src.write_text(imports_only(text))
    r = subprocess.run([sys.executable, "-c", RUNNER, str(src), str(fake)], cwd=ROOT, capture_output=True, text=True,
                       timeout=300, env={"PATH": "/usr/bin:/bin", "MPLBACKEND": "Agg", "OMP_NUM_THREADS": "1",
                                         "HOME": str(Path.home())})
    err = r.stderr.strip().splitlines()[-1] if r.returncode else ""
    return r.returncode == 0, err


def check(item, rev, tmp):
    name, d = item
    old = subprocess.run(["git", "show", f"{rev}:scripts/{name}.py"], capture_output=True, text=True, check=True).stdout
    new = (ROOT / "scripts" / d / f"{name}.py").read_text()
    ok_old, e_old = run(old, ROOT / "scripts" / f"{name}.py", "old", name, tmp)
    ok_new, e_new = run(new, ROOT / "scripts" / d / f"{name}.py", "new", name, tmp)
    return name, d, ok_old, ok_new, e_old, e_new


def check_in_place(path, rev, tmp):
    """A file that did not move but whose imports were rewritten (--also): old text vs new text, same place."""
    p = Path(path)
    old = subprocess.run(["git", "show", f"{rev}:{path}"], capture_output=True, text=True, check=True).stdout
    tag = p.with_suffix("").as_posix().replace("/", "__")
    ok_old, e_old = run(old, ROOT / p, "old", tag, tmp)
    ok_new, e_new = run((ROOT / p).read_text(), ROOT / p, "new", tag, tmp)
    return p.stem, str(p.parent), ok_old, ok_new, e_old, e_new


def main():
    moves = json.loads(Path(sys.argv[1]).read_text())
    rev = sys.argv[2]
    also = sys.argv[sys.argv.index("--also") + 1:] if "--also" in sys.argv else []
    argv = sys.argv[:sys.argv.index("--also")] if also else sys.argv
    workers = int(argv[3]) if len(argv) > 3 else 4
    tmp = Path(subprocess.run(["mktemp", "-d"], capture_output=True, text=True).stdout.strip())
    with ThreadPoolExecutor(workers) as ex:
        res = list(ex.map(lambda it: check(it, rev, tmp), sorted(moves.items())))
        res += list(ex.map(lambda f: check_in_place(f, rev, tmp), also))
    reg = [r for r in res if r[2] and not r[3]]
    both_bad = [r for r in res if not r[2] and not r[3]]
    fixed = [r for r in res if not r[2] and r[3]]
    print(f"checked {len(res)}: ok at both {sum(1 for r in res if r[2] and r[3])}, REGRESSIONS {len(reg)}, "
          f"failing at both (environment, not the move) {len(both_bad)}, failing only at the old place {len(fixed)}")
    for r in reg:
        print(f"  REGRESSION {r[1]}/{r[0]}: {r[5]}")
    for r in both_bad:
        print(f"  both {r[1]}/{r[0]}: old: {r[4]} | new: {r[5]}")
    sys.exit(1 if reg else 0)


if __name__ == "__main__":
    main()
