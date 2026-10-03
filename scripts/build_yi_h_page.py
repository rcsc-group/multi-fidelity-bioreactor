"""Build the Yi-h artifact page: the guide, then the specification in a closed section.

Renders docs/yi_h_guide.md + docs/mf_problem_definition.md to static HTML on the
server (python-markdown), so the viewer's native Mermaid renderer and MathJax both
see static content. Math is protected from the markdown pass and restored verbatim.
Figures are referenced as fig/<path under experiments/multifidelity>.
Published at https://claude.ai/artifact/9i8VX83AGq7CLm5wFNqdzi
"""
import html
import re
from pathlib import Path

import markdown

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "build/yi_h_page/yi_h.html"


def load():
    guide = (REPO / "docs/yi_h_guide.md").read_text()
    spec = (REPO / "docs/mf_problem_definition.md").read_text()
    spec = re.sub(r"^(#{2,5}) ", lambda m: "#" + m.group(1) + " ", spec, flags=re.M)   # demote
    spec = re.sub(r"^# .*$", "## Specification", spec, count=1, flags=re.M)
    head, body = spec.split("## Specification", 1)
    spec = ("## Specification\n\n<details markdown=\"1\">\n<summary>The precise definition, "
            "algorithm and assumptions (for implementation and review)</summary>\n\n" + body
            + "\n\n</details>\n")
    text = guide + "\n\n---\n\n" + spec
    return re.sub(r"\]\(\.\./experiments/multifidelity/([^)]+)\)", r"](fig/\1)", text)


def render(text):
    math = []

    def stash(m):
        math.append(m.group(0))
        return f"MATHTOKEN{len(math) - 1}X"
    text = re.sub(r"\$\$[\s\S]+?\$\$", stash, text)
    text = re.sub(r"\\\([\s\S]+?\\\)", stash, text)
    out = markdown.markdown(text, extensions=["tables", "fenced_code", "md_in_html"])
    out = re.sub(r"MATHTOKEN(\d+)X", lambda m: html.escape(math[int(m.group(1))], quote=False), out)
    out = re.sub(r'<pre><code class="language-mermaid">([\s\S]*?)</code></pre>',
                 lambda m: '<pre class="mermaid">' + html.unescape(m.group(1)) + "</pre>", out)
    out = re.sub(r"\[(lit|data|inference|assumption)([^\]<]*)\]",
                 lambda m: f'<span class="chip chip-{m.group(1)}">{m.group(1)}{m.group(2)}</span>', out)
    out = re.sub(r"<table>", '<div class="table-wrap"><table>', out)
    out = re.sub(r"</table>", "</table></div>", out)
    out = re.sub(r"<p>(<img [^>]+>)</p>", r'<figure class="fig">\1</figure>', out)
    toc, n = [], [0]

    def h2(m):
        n[0] += 1
        title = re.sub(r"<[^>]+>", "", m.group(1))
        toc.append(f'<li><a href="#s{n[0]}">{title}</a></li>')
        return f'<h2 id="s{n[0]}">{m.group(1)}</h2>'
    out = re.sub(r"<h2>(.*?)</h2>", h2, out)
    return out, "\n".join(toc)


STYLE = (REPO / "scripts/yi_h_page_style.css").read_text() if (REPO / "scripts/yi_h_page_style.css").exists() else ""


def main():
    body, toc = render(load())
    page = f"""<title>Yi-h Problem Statement</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans+Condensed:wght@500;600;700&family=Source+Serif+4:ital,opsz,wght@0,8..60,400;0,8..60,600;1,8..60,400&display=swap">
<style>
{STYLE}
</style>
<div class="shell">
  <nav class="toc" aria-label="Sections"><p>Sections</p><ol>
{toc}
  </ol></nav>
  <main id="doc">
{body}
  </main>
</div>
<script>
window.MathJax = {{ tex: {{ inlineMath: [["\\\\(", "\\\\)"]], displayMath: [["$$", "$$"]] }}, svg: {{ fontCache: "global" }} }};
</script>
<script src="https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-svg.js"></script>
"""
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(page)
    print("wrote", OUT, OUT.stat().st_size, "bytes;", len(re.findall("<h2", page)), "sections")


if __name__ == "__main__":
    main()
