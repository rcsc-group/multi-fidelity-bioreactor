"""Build the Yi-h artifact page from docs/mf_problem_definition.md."""
import html
import re
from pathlib import Path

REPO = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
OUT = REPO / "build/yi_h_page/yi_h.html"  # artifact https://claude.ai/artifact/9i8VX83AGq7CLm5wFNqdzi; publish with fig/*.png from experiments/multifidelity

md = (REPO / "docs/mf_problem_definition.md").read_text()
md = re.sub(r"\]\(\.\./experiments/multifidelity/([^)]+)\)", r"](fig/\1)", md)
assert "</script" not in md

PAGE = r"""<title>Yi-h Problem Statement</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans+Condensed:wght@500;600;700&family=Source+Serif+4:ital,opsz,wght@0,8..60,400;0,8..60,600;1,8..60,400&display=swap">
<style>
/* Layout: one reading column (~70ch) with a sticky section index on wide screens; drawing-title-block headings. */
:root {
  --ground: #f5f7f6;
  --paper: #ffffff;
  --ink: #18221f;
  --ink-soft: #4a5753;
  --rule: #d5ddd9;
  --accent: #2b5a87;
  --chip-lit: #2e7056;
  --chip-data: #2b5a87;
  --chip-inf: #6a4d9c;
  --chip-ass: #9a6512;
  --code-bg: #eaefed;
  --font-display: "IBM Plex Sans Condensed", "Arial Narrow", "Helvetica Neue", Arial, sans-serif;
  --font-body: "Source Serif 4", Georgia, "Times New Roman", serif;
  --font-mono: "IBM Plex Mono", ui-monospace, Menlo, Consolas, monospace;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --ground: #101614; --paper: #161e1b; --ink: #dce5e1; --ink-soft: #9fb0aa; --rule: #2c3833;
    --accent: #86b2dc; --chip-lit: #7cc4a4; --chip-data: #86b2dc; --chip-inf: #b9a0e3; --chip-ass: #e0ae5c;
    --code-bg: #1f2925; color-scheme: dark;
  }
}
:root[data-theme="dark"] {
  --ground: #101614; --paper: #161e1b; --ink: #dce5e1; --ink-soft: #9fb0aa; --rule: #2c3833;
  --accent: #86b2dc; --chip-lit: #7cc4a4; --chip-data: #86b2dc; --chip-inf: #b9a0e3; --chip-ass: #e0ae5c;
  --code-bg: #1f2925; color-scheme: dark;
}
body { background: var(--ground); color: var(--ink); font-family: var(--font-body); font-size: 17px; line-height: 1.6; }
.shell { display: grid; grid-template-columns: minmax(0, 1fr); gap: 2.5rem; max-width: 1180px; margin: 0 auto; padding-inline: 16px; padding-block: 2rem 4rem; }
@media (min-width: 1000px) { .shell { grid-template-columns: 220px minmax(0, 1fr); padding-inline: 32px; } }
nav.toc { display: none; }
@media (min-width: 1000px) {
  nav.toc { display: block; position: sticky; top: calc(env(safe-area-inset-top, 0px) + 1.5rem); align-self: start; max-height: calc(100vh - 3rem); overflow-y: auto; }
}
nav.toc p { font-family: var(--font-display); font-size: 0.75rem; letter-spacing: 0.08em; text-transform: uppercase; color: var(--ink-soft); margin: 0 0 0.6rem; }
nav.toc ol { list-style: none; margin: 0; padding: 0; display: grid; gap: 0.35rem; border-left: 2px solid var(--rule); }
nav.toc a { display: block; padding-left: 0.75rem; font-family: var(--font-display); font-size: 0.92rem; color: var(--ink-soft); text-decoration: none; line-height: 1.3; }
nav.toc a:hover, nav.toc a:focus-visible { color: var(--accent); }
main { min-width: 0; max-width: 46rem; }
main > * { min-width: 0; }
h1, h2, h3 { font-family: var(--font-display); line-height: 1.2; text-wrap: balance; color: var(--ink); }
h1 { font-size: 2.1rem; font-weight: 700; margin: 0 0 0.75rem; letter-spacing: -0.005em; }
h2 { font-size: 1.5rem; font-weight: 600; margin: 3rem 0 0.9rem; padding-top: 0.6rem; border-top: 2px solid var(--ink); }
h3 { font-size: 1.15rem; font-weight: 600; margin: 2rem 0 0.6rem; color: var(--accent); }
p, li { max-width: 70ch; }
a { color: var(--accent); }
hr { border: 0; border-top: 1px solid var(--rule); margin: 2.5rem 0; }
code { font-family: var(--font-mono); font-size: 0.85em; background: var(--code-bg); padding: 0.08em 0.3em; border-radius: 3px; overflow-wrap: anywhere; }
strong { font-weight: 600; }
ul, ol { padding-left: 1.4rem; }
li + li { margin-top: 0.25rem; }
.table-wrap { overflow-x: auto; margin: 1.2rem 0; border: 1px solid var(--rule); background: var(--paper); }
table { border-collapse: collapse; width: 100%; font-size: 0.9rem; line-height: 1.45; }
th, td { text-align: left; vertical-align: top; padding: 0.55rem 0.7rem; border-bottom: 1px solid var(--rule); }
th { font-family: var(--font-display); font-weight: 600; font-size: 0.85rem; letter-spacing: 0.02em; background: var(--code-bg); }
tr:last-child td { border-bottom: 0; }
td:first-child { white-space: nowrap; }
figure.fig { margin: 1.6rem 0 0.4rem; background: var(--paper); border: 1px solid var(--rule); padding: 0.75rem; }
figure.fig img { display: block; width: 100%; height: auto; background: #ffffff; }
p em:only-child { color: var(--ink-soft); font-size: 0.92rem; }
.chip { display: inline-block; font-family: var(--font-mono); font-size: 0.72em; font-weight: 500; line-height: 1; padding: 0.22em 0.45em 0.2em; border-radius: 2px; border: 1px solid currentColor; vertical-align: 0.12em; white-space: nowrap; }
.chip-lit { color: var(--chip-lit); }
.chip-data { color: var(--chip-data); }
.chip-inference { color: var(--chip-inf); }
.chip-assumption { color: var(--chip-ass); }
mjx-container[display="true"] { overflow-x: auto; overflow-y: hidden; padding: 0.3rem 0; }
.legend { display: flex; flex-wrap: wrap; gap: 0.5rem 1rem; font-size: 0.85rem; color: var(--ink-soft); margin: 0.4rem 0 1.5rem; }
@media (prefers-reduced-motion: reduce) { * { scroll-behavior: auto !important; } }
</style>

<div class="shell">
  <nav class="toc" aria-label="Sections"><p>Sections</p><ol id="toc"></ol></nav>
  <main id="doc"><p>Loading the document…</p></main>
</div>

<script type="text/markdown" id="src">
__MD__
</script>
<script src="https://cdn.jsdelivr.net/npm/marked@12.0.2/marked.min.js"></script>
<script>
window.MathJax = {
  tex: { inlineMath: [["\\(", "\\)"]], displayMath: [["$$", "$$"]] },
  svg: { fontCache: "global" },
  startup: { typeset: false }
};
</script>
<script src="https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-svg.js"></script>
<script>
(function () {
  var src = document.getElementById("src").textContent;
  var math = [];
  function stash(m) { math.push(m); return "@@M" + (math.length - 1) + "@@"; }
  src = src.replace(/\$\$[\s\S]+?\$\$/g, stash).replace(/\\\([\s\S]+?\\\)/g, stash);
  var out = marked.parse(src, { gfm: true });
  var esc = function (s) { return s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;"); };
  out = out.replace(/@@M(\d+)@@/g, function (_, i) { return esc(math[+i]); });
  out = out.replace(/\[(lit|data|inference|assumption)([^\]<]*)\]/g, function (_, kind, rest) {
    return '<span class="chip chip-' + kind + '">' + kind + rest + "</span>";
  });
  var doc = document.getElementById("doc");
  doc.innerHTML = out;
  doc.querySelectorAll("table").forEach(function (t) {
    var w = document.createElement("div"); w.className = "table-wrap";
    t.parentNode.insertBefore(w, t); w.appendChild(t);
  });
  doc.querySelectorAll("p > img:only-child").forEach(function (img) {
    var p = img.parentNode, f = document.createElement("figure"); f.className = "fig";
    img.loading = "lazy"; p.parentNode.replaceChild(f, p); f.appendChild(img);
  });
  var toc = document.getElementById("toc");
  doc.querySelectorAll("h2").forEach(function (h, i) {
    h.id = "s" + i;
    var li = document.createElement("li"), a = document.createElement("a");
    a.href = "#s" + i; a.textContent = h.textContent; li.appendChild(a); toc.appendChild(li);
  });
  function typeset() { if (window.MathJax && MathJax.typesetPromise) { MathJax.typesetPromise([doc]); } }
  if (window.MathJax && MathJax.startup && MathJax.startup.promise) { MathJax.startup.promise.then(typeset); }
  else { window.addEventListener("load", typeset); }
})();
</script>
"""

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(PAGE.replace("__MD__", md))
print("wrote", OUT, OUT.stat().st_size, "bytes")
