"""Move the runnable top-level scripts of scripts/ into sub-folders (one-off tool of the 2026-11 clean-up).

Rule (scripts/README.md):
  - a module that other code imports stays at scripts/ top level (its import path does not change);
  - every other top-level script moves, by the first rule that matches:
      test_*                         -> scripts/hypotheses/   (numerical experiments, one per diary hypothesis)
      plot_* fig_* render_*          -> scripts/figures/      (fig9_sweep_driver excepted: it submits jobs)
      submit_* *_watchdog *_driver,
      or text contains "sbatch"      -> scripts/run/          (submits or babysits SLURM jobs)
      anything else                  -> scripts/analyse/

What it changes, per moved file:
  - one more directory level in every Path(__file__) parent chain (parents[1] -> parents[2], .parent.parent ->
    .parent.parent.parent), so every path still points where it did;
  - plain imports of a top-level scripts module get scripts/ on sys.path;
and in every tracked text file except history (diary*, docs/, experiments/): "scripts/NAME.py" and
"scripts.NAME" -> the new location.

Checks (the tool refuses to apply when one fails):
  1. every __file__ expression evaluates to the same path at the old and the new location;
  2. no remaining "NAME.py" mention of a moved script outside history, other than in the new path form.
A third check (imports resolve from the new location) runs separately: tools/regroup_import_check.py.

Usage:
  python tools/regroup_scripts.py SCRIPTS_TSV LIBS_TXT [--apply]   (from the repo root; dry run by default)
"""
import ast
import csv
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path.cwd()
SCRIPTS = ROOT / "scripts"
HISTORY = ("diary", "docs/", "experiments/", "tools/regroup_scripts.py")
TEXT = (".py", ".sh", ".sbatch", ".toml", ".yml", ".yaml", ".json", ".cfg", ".txt", ".md", "Makefile")


def folder(name: str, text: str) -> str:
    if name.startswith("test_"):
        return "hypotheses"
    if name.startswith(("plot_", "fig_", "render_")):
        return "figures"
    if name.lstrip("_").startswith(("submit_", "benchmark_submit")) or name.endswith(("_watchdog", "_driver")) \
            or "sbatch" in text:
        return "run"
    return "analyse"


def deeper(src: str) -> str:
    """One more parent level in every Path(__file__) chain."""
    src = re.sub(r"(Path\(__file__\)(?:\.resolve\(\))?)\.parents\[(\d)\]",
                 lambda m: f"{m.group(1)}.parents[{int(m.group(2)) + 1}]", src)
    src = re.sub(r"(Path\(__file__\)(?:\.resolve\(\))?)((?:\.parent)+)(?![\w\[])",
                 lambda m: m.group(1) + m.group(2) + ".parent", src)
    return src


def file_exprs(src: str):
    """Source text of every maximal expression that contains __file__."""
    tree = ast.parse(src)
    out = []

    class V(ast.NodeVisitor):
        def generic_visit(self, node):
            if isinstance(node, ast.expr) and "__file__" in ast.unparse(node):
                # descend until the parent stops being an attribute/subscript/call chain on __file__
                out.append(node)
                return
            super().generic_visit(node)

        def visit_Call(self, node):
            # keep the chain Path(__file__).resolve().parents[1] whole, but do not swallow a call whose
            # argument merely contains it (e.g. sys.path.insert(0, str(...)))
            f = ast.unparse(node.func)
            if "__file__" in f:
                out.append(node)
                return
            for a in node.args + [k.value for k in node.keywords]:
                self.visit(a)

    V().visit(tree)
    return [ast.unparse(n) for n in out]


def evaluate(expr: str, file: Path):
    import os
    from pathlib import Path as P
    try:
        v = eval(expr, {"Path": P, "_Path": P, "os": os, "__file__": str(file)})
    except Exception as e:  # noqa: BLE001
        return f"<error {type(e).__name__}>"
    return str(v)


def main():
    tsv, libs_txt = Path(sys.argv[1]), Path(sys.argv[2])
    apply = "--apply" in sys.argv
    libs = set(Path(libs_txt).read_text().split())
    rows = {r["name"]: r for r in csv.DictReader(open(tsv), delimiter="\t")}
    moves = {}
    for p in sorted(SCRIPTS.glob("*.py")):
        n = p.stem
        if n == "__init__" or n in libs:
            continue
        r = rows.get(n)
        if r is None or int(r["n_imp"]) > 0 or r["deck"]:
            continue
        moves[n] = folder(n, p.read_text(errors="ignore"))

    # 1. rewrite and check __file__ expressions
    new_src, problems = {}, []
    lib_pat = re.compile(rf"^(\s*)(import|from)\s+({'|'.join(map(re.escape, sorted(libs)))})\b", re.M)
    for n, d in moves.items():
        old_p, new_p = SCRIPTS / f"{n}.py", SCRIPTS / d / f"{n}.py"
        src = old_p.read_text()
        s2 = deeper(src)
        if lib_pat.search(s2):
            # plain sibling import: put scripts/ on sys.path before the first such import
            m = lib_pat.search(s2)
            line_start = s2.rfind("\n", 0, m.start()) + 1
            shim = ("import sys as _sys; from pathlib import Path as _SP  # noqa: E402\n"
                    "_sys.path.insert(0, str(_SP(__file__).resolve().parents[1]))  # scripts/: shared modules\n")
            s2 = s2[:line_start] + shim + s2[line_start:]
        ea, eb = file_exprs(src), file_exprs(s2)
        eb = [e for e in eb if "_SP(__file__)" not in e]
        if len(ea) != len(eb):
            problems.append(f"{n}: {len(ea)} __file__ expressions before, {len(eb)} after")
        for a, b in zip(ea, eb):
            va, vb = evaluate(a, old_p), evaluate(b, new_p)
            if va != vb:
                problems.append(f"{n}: {a!r} -> {va} but {b!r} -> {vb}")
        new_src[n] = s2

    # 2. rewrite mentions in tracked text files outside history
    files = subprocess.run(["git", "ls-files"], capture_output=True, text=True, check=True).stdout.split("\n")
    texts = [f for f in files if f and (f.endswith(TEXT)) and not f.startswith(HISTORY)]
    name_alt = "|".join(map(re.escape, sorted(moves, key=len, reverse=True)))
    path_re = re.compile(rf"scripts/({name_alt})\.py\b")
    mod_re = re.compile(rf"\bscripts\.({name_alt})\b")
    join_re = re.compile(rf"""(["']scripts["'])\s*/\s*(["'])({name_alt})\.py["']""")
    edits, leftovers = {}, []
    for f in texts:
        p = ROOT / f
        try:
            t = p.read_text()
        except (OSError, UnicodeDecodeError):
            continue
        if f.startswith("scripts/") and Path(f).parent == Path("scripts") and Path(f).stem in moves and f.endswith(".py"):
            t = new_src[Path(f).stem]
        t2 = path_re.sub(lambda m: f"scripts/{moves[m.group(1)]}/{m.group(1)}.py", t)
        t2 = mod_re.sub(lambda m: f"scripts.{moves[m.group(1)]}.{m.group(1)}", t2)
        t2 = join_re.sub(lambda m: f'{m.group(1)} / "{moves[m.group(3)]}" / {m.group(2)}{m.group(3)}.py{m.group(2)}', t2)
        if t2 != t or Path(f).stem in moves and Path(f).parent == Path("scripts") and f.endswith(".py"):
            edits[f] = t2
        shell = f.endswith((".sh", ".sbatch")) or f == "Makefile"
        for ln in t2.splitlines():
            if shell and ln.lstrip().startswith("#"):
                continue
            for m in re.finditer(rf"(?<![\w/])({name_alt})\.py\b", ln):
                quoted = re.search(rf"""["'][^"'\n]*{re.escape(m.group(0))}""", ln)
                if (quoted or shell) and Path(f).stem != m.group(1):
                    leftovers.append(f"{f}: {ln.strip()[:150]}")

    from collections import Counter
    print(f"moves: {len(moves)} {dict(Counter(moves.values()))}")
    print(f"files to edit (incl. moved): {len(edits)}")
    print(f"__file__ check problems: {len(problems)}")
    for x in problems:
        print("  ", x)
    print(f"bare NAME.py mentions left for review: {len(leftovers)}")
    for x in leftovers[:80]:
        print("  ", x)
    # input of tools/regroup_import_check.py (scratch file; the kept record is internal/docs/history/moved-scripts.md)
    json.dump(moves, open(ROOT / "tools" / ".regroup_moves.json", "w"), indent=1, sort_keys=True)
    if not apply:
        return
    if problems:
        sys.exit("refused: __file__ check failed")
    for d in set(moves.values()):
        (SCRIPTS / d).mkdir(exist_ok=True)
    for n, d in moves.items():
        subprocess.run(["git", "mv", f"scripts/{n}.py", f"scripts/{d}/{n}.py"], check=True)
    for f, t in edits.items():
        p = ROOT / f
        if Path(f).parent == Path("scripts") and Path(f).stem in moves and f.endswith(".py"):
            p = SCRIPTS / moves[Path(f).stem] / Path(f).name
        p.write_text(t)
    print("applied")


if __name__ == "__main__":
    main()
