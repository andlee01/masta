# Editing this chapter

Read [Equation Generation](equation_generation.md) on GitHub. Edit
[equation_generation.qmd](equation_generation.qmd) and regenerate the Markdown
from the repository root:

```bash
quarto render doc/equation_generation/equation_generation.qmd --to gfm
```

Commit both the `.qmd` source and generated `.md` file. The source also supports
HTML output with `--to html`.

The chapter was migrated from `~/website/html/equation_generation.html` on
2026-11-06. Its 28 figure assets are SVGs in `figures/`; the corresponding 25
Graphviz `.dot` files and three standalone Circuitikz `.tex` files are in
`sources/`. The legacy originals were preserved. Windows download metadata,
website menus, CSS, and MathJax scripts were excluded.

Rebuild all figures from the repository root:

```bash
python3 scripts/generate_figures_equation_generation.py
quarto render doc/equation_generation/equation_generation.qmd --to gfm
```

The renderer requires Graphviz (`dot`), LaTeX with Circuitikz and its dependencies,
and `dvisvgm`. It uses `latex` to generate DVI and `dvisvgm --no-fonts` to generate
vector paths, including text, with no embedded raster circuit image. Intermediate
files live in a temporary directory. Use `--tex-only` to regenerate just the
three circuit diagrams if Graphviz is unavailable.

The legacy script was `~/website/scripts/generate_figures_equation_generation.sh`.
It used `dot -Tsvg` for graphs and `pdflatex` followed by ImageMagick `convert`
for circuit images. The new renderer uses the chapter-local sources; `fig1.dot`
still maps to `fig1a.svg`, and circuit TeX files map to `circuitN_tex.svg`.

Commit the sources, generated SVGs, and generated Markdown together. GitHub's
repository Markdown view displays the referenced images but does not execute
Graphviz, LaTeX, or Quarto. Reading the chapter requires none of those tools.
PDF, DVI, auxiliary files, and temporary PNG previews need not be committed.

Multi-panel figures use HTML table rows that work in GitHub Markdown. Circuit
panels are 480 pixels wide; graph panels are 200 pixels high with proportional
widths, so tall, narrow graphs do not overwhelm the circuit drawings.

On the current Linux/WSL computer, Graphviz is installed per user under
`~/.local/opt/graphviz`, with a `~/.local/bin/dot` launcher. This avoids requiring
a sudo password. Other computers can use their normal Graphviz installation;
the renderer only requires `dot` to be on `PATH`.

Migration retains the legacy exposition and equations rather than updating
them to the current Python API. Equation labels were replaced with explicit
Markdown links, figure groups were numbered, and minor typographical errors
were corrected. The basic-loop definition was corrected to count co-tree
branches. The remaining derivations have not received a full technical review.
