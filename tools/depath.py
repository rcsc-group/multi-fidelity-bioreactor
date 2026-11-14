"""Replace hard-coded OSCAR folders in Python files with scripts/paths.py (one-off tool of the 2026-11 clean-up).

    python tools/depath.py FILE ... [--apply]        (from the repo root; dry run by default)

Python files are rewritten as below; .sh/.sbatch files by rewrite_shell().

Per file:
  - a string literal that holds a folder becomes an f-string: "/oscar/scratch/eaguerov/x" -> f"{paths.SCRATCH}/x";
    other braces are doubled. CHECK: with the default (OSCAR) folders, the new literal evaluates to the old string,
    character for character, or the file is refused.
  - an f-string: the folder text in it becomes {paths.X}.
  - a docstring or a comment: the folder text becomes $MFB_SCRATCH etc. (for the reader; no effect on the code).
  - `from scripts import paths` is added before the first use: after an existing `from scripts... import`, or
    else with a repo-root sys.path line (the file may run as `python scripts/x/y.py`).
A bytes literal, or a folder inside an f-string expression, refuses the file (fix by hand).
"""
import ast
import io
import re
import sys
import tokenize
from pathlib import Path

ROOT = Path.cwd()
# longest first; (folder, attribute, variable for human text)
MAP = [("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor", "REPO", "$MFB_REPO"),
       ("/oscar/data/dharri15/eaguerov/Github/BioReactor3D/dev/rocking-bioreactor-2d", "LEGACY_2D", "$MFB_LEGACY_2D"),
       ("/oscar/data/dharri15/eaguerov/basilisk/src", "BASILISK", "$BASILISK"),
       ("/oscar/data/dharri15/eaguerov", "DATA", "$MFB_DATA"),
       ("/oscar/scratch/eaguerov", "SCRATCH", "$MFB_SCRATCH")]
DEFAULTS = {a: f for f, a, _ in MAP}
# (?<!:-): a folder that is already the default of ${VAR:-...} is left alone, so a second pass changes nothing
PAT = re.compile(r"(?<!:-)(?:" + "|".join(re.escape(f) for f, _, _ in MAP) + r")(?![\w-])")
ATTR = {f: a for f, a, _ in MAP}
HUMAN = {f: h for f, _, h in MAP}


class _P:  # the default folders, as scripts.paths gives them with no environment set
    pass


for a, f in DEFAULTS.items():
    setattr(_P, a, Path(f))


def docstring_starts(tree) -> set:
    out = set()
    for n in ast.walk(tree):
        if isinstance(n, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) and n.body:
            s = n.body[0]
            if isinstance(s, ast.Expr) and isinstance(s.value, ast.Constant) and isinstance(s.value.value, str):
                out.add((s.lineno, s.col_offset))
    return out


def rewrite(src: str, rel: Path):
    lines = src.splitlines(keepends=True)
    off = [0]
    for ln in lines:
        off.append(off[-1] + len(ln))
    pos = lambda lc: off[lc[0] - 1] + lc[1]  # noqa: E731
    tree = ast.parse(src)
    docs = docstring_starts(tree)
    edits, problems, first_code_use = [], [], None
    toks = list(tokenize.generate_tokens(io.StringIO(src).readline))
    i = 0
    while i < len(toks):
        t = toks[i]
        if t.type == tokenize.COMMENT and PAT.search(t.string):
            edits.append((pos(t.start), pos(t.end), PAT.sub(lambda m: HUMAN[m.group(0)], t.string)))
        elif t.type == tokenize.STRING and PAT.search(t.string):
            if t.start in docs:
                edits.append((pos(t.start), pos(t.end), PAT.sub(lambda m: HUMAN[m.group(0)], t.string)))
            else:
                m = re.match(r"([A-Za-z]*)(.*)", t.string, re.S)
                pre, body = m.group(1), m.group(2)
                if "b" in pre.lower():
                    problems.append(f"bytes literal at line {t.start[0]}")
                else:
                    q = 3 if body[:3] in ('"""', "'''") else 1
                    inner = body[q:-q].replace("{", "{{").replace("}", "}}")
                    inner = PAT.sub(lambda mm: "{paths." + ATTR[mm.group(0)] + "}", inner)
                    new = pre.replace("u", "").replace("U", "") + "f" + body[:q] + inner + body[-q:]
                    try:
                        same = eval(new, {"paths": _P}) == ast.literal_eval(t.string)
                    except Exception as e:  # noqa: BLE001
                        same, new = False, f"{new}  [{e}]"
                    if not same:
                        problems.append(f"line {t.start[0]}: {t.string[:60]} -> {new[:80]} is not the same string")
                    edits.append((pos(t.start), pos(t.end), new))
                    first_code_use = min(first_code_use or 10**9, t.start[0])
        elif t.type == tokenize.FSTRING_START:
            j, depth = i + 1, 1
            while depth:
                depth += toks[j].type == tokenize.FSTRING_START
                depth -= toks[j].type == tokenize.FSTRING_END
                j += 1
            a, b = pos(t.start), pos(toks[j - 1].end)
            seg = src[a:b]
            if PAT.search(seg):
                mids = [k for k in toks[i:j] if k.type == tokenize.FSTRING_MIDDLE]
                if not any(PAT.search(k.string) for k in mids):
                    problems.append(f"folder inside an f-string expression at line {t.start[0]}")
                edits.append((a, b, PAT.sub(lambda mm: "{paths." + ATTR[mm.group(0)] + "}", seg)))
                first_code_use = min(first_code_use or 10**9, t.start[0])
            i = j
            continue
        i += 1
    if not edits:
        return None, problems
    new = src
    for a, b, s in sorted(edits, reverse=True):
        new = new[:a] + s + new[b:]
    if first_code_use is not None and not re.search(r"^from scripts import paths\b", new, re.M):
        new = add_import(new, tree, first_code_use, rel)
    return new, problems


def add_import(new: str, tree, first_use: int, rel: Path) -> str:
    """Add `from scripts import paths` on its own line before line `first_use` (a top-level statement boundary)."""
    lines = new.splitlines(keepends=True)
    tops = [s for s in tree.body if s.end_lineno < first_use]
    from_scripts = [s for s in tops if isinstance(s, ast.ImportFrom) and (s.module or "").split(".")[0] == "scripts"]
    if from_scripts:
        at = from_scripts[0].end_lineno
        block = "from scripts import paths  # noqa: E402  (machine-specific folders)\n"
    else:
        imps = [s for s in tops if isinstance(s, (ast.Import, ast.ImportFrom))]
        doc = tree.body[0] if tree.body and isinstance(tree.body[0], ast.Expr) and \
            isinstance(getattr(tree.body[0], "value", None), ast.Constant) else None
        k = len(rel.parts) - 1
        # a top-level sys.path line before the first use that puts the repo root first: the import alone is enough
        root_path = [s for s in tops if "sys.path" in ast.unparse(s) and
                     (f"parents[{k}]" in ast.unparse(s) or "ROOT" in ast.unparse(s))]
        names = {a.asname or a.name for s in imps for a in s.names}
        if root_path:
            at = max(s.end_lineno for s in root_path + imps)
            block = "from scripts import paths  # noqa: E402  (machine-specific folders)\n"
        else:
            at = imps[-1].end_lineno if imps else (doc.end_lineno if doc else 0)
            block = ("" if "sys" in names else "import sys  # noqa: E402\n") + \
                    ("" if "Path" in names else "from pathlib import Path  # noqa: E402\n") + \
                    f"sys.path.insert(0, str(Path(__file__).resolve().parents[{k}]))  # repo root\n" + \
                    "from scripts import paths  # noqa: E402  (machine-specific folders)\n"
    return "".join(lines[:at]) + block + "".join(lines[at:])


def rewrite_shell(src: str):
    """Shell/SLURM: code -> ${MFB_X:-<OSCAR default>}; comments -> $MFB_X; #SBATCH lines unchanged (SLURM does
    not expand variables there). A folder inside single quotes refuses the file (no expansion there)."""
    out, problems = [], []
    for n, ln in enumerate(src.splitlines(keepends=True), 1):
        s = ln.lstrip()
        if s.startswith("#SBATCH") or not PAT.search(ln):
            out.append(ln)
        elif s.startswith("#"):
            out.append(PAT.sub(lambda m: HUMAN[m.group(0)], ln))
        else:
            if any(PAT.search(q) for q in re.findall(r"'[^']*'", ln)):
                problems.append(f"line {n}: folder inside single quotes")
            out.append(PAT.sub(lambda m: "${" + HUMAN[m.group(0)][1:] + ":-" + m.group(0) + "}", ln))
    new = "".join(out)
    return (new if new != src else None), problems


def main():
    files = [Path(a) for a in sys.argv[1:] if not a.startswith("--")]
    apply = "--apply" in sys.argv
    bad, done = {}, []
    for f in files:
        shell = f.suffix in (".sh", ".sbatch")
        new, problems = (rewrite_shell if shell else rewrite)(*((f.read_text(),) if shell else (f.read_text(), f)))
        if problems:
            bad[str(f)] = problems
        elif new is not None:
            done.append((f, new))
    print(f"files: {len(files)}  to rewrite: {len(done)}  refused: {len(bad)}")
    for f, p in bad.items():
        print(f"  REFUSED {f}: {p[:3]}")
    if apply:
        for f, new in done:
            f.write_text(new)
        print("applied")


if __name__ == "__main__":
    main()
