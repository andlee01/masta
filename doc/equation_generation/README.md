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
2026-11-06. All 28 unique referenced images are in `figures/`; the corresponding
25 Graphviz `.dot` files and three standalone Circuitikz `.tex` files are in
`sources/`. The legacy originals were preserved. Windows download metadata,
website menus, CSS, and MathJax scripts were excluded.

The Graphviz source `fig1.dot` produces `figures/fig1a.svg`; other `.dot` files
match their SVG stems. For example:

```bash
dot -Tsvg doc/equation_generation/sources/fig1.dot -o doc/equation_generation/figures/fig1a.svg
```

The circuit PNGs came from the matching `circuit1.tex`, `circuit2.tex`, and
`circuit3.tex` sources. Rebuilding those requires LaTeX with Circuitikz and a
PDF-to-image converter; reading the chapter requires neither.

Migration retains the legacy exposition and equations rather than updating
them to the current Python API. Equation labels were replaced with explicit
Markdown links, figure groups were numbered, and minor typographical errors
were corrected. The basic-loop definition was corrected to count co-tree
branches. The remaining derivations have not received a full technical review.
