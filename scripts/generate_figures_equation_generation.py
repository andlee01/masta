"""Render the equation-generation drawings as vector SVGs.

Run from any directory. Requires Graphviz (dot), LaTeX with Circuitikz, and
dvisvgm. Use --tex-only to rebuild just the three circuit drawings.
"""

import argparse
from pathlib import Path
import shutil
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tex-only", action="store_true")
    args = parser.parse_args()
    chapter = Path(__file__).resolve().parents[1] / "doc/equation_generation"
    sources = chapter / "sources"
    figures = chapter / "figures"
    required = ["latex", "dvisvgm"] + ([] if args.tex_only else ["dot"])
    missing = [tool for tool in required if shutil.which(tool) is None]
    if missing:
        parser.error("Missing required tools: " + ", ".join(missing))

    # Build everything in isolation before replacing any published figures.
    with tempfile.TemporaryDirectory(prefix="masta-equation-figures-") as tmp:
        work = Path(tmp)
        outputs = []
        for source in sorted(sources.glob("*.tex")):
            result = subprocess.run(
                ["latex", "-interaction=nonstopmode", "-halt-on-error",
                 "-no-shell-escape", f"-output-directory={work}", str(source)],
                cwd=work, text=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
            if result.returncode:
                raise RuntimeError(f"LaTeX failed for {source.name}:\n{result.stdout}")
            output = work / f"{source.stem}_tex.svg"
            subprocess.run(
                ["dvisvgm", "--no-fonts", "--exact-bbox",
                 f"--output={output}", str(work / f"{source.stem}.dvi")],
                cwd=work, check=True,
            )
            outputs.append(output)
        if not args.tex_only:
            for source in sorted(sources.glob("*.dot")):
                stem = "fig1a" if source.stem == "fig1" else source.stem
                output = work / f"{stem}.svg"
                subprocess.run(["dot", "-Tsvg", str(source), "-o", str(output)],
                               check=True)
                outputs.append(output)
        figures.mkdir(exist_ok=True)
        for output in outputs:
            shutil.copyfile(output, figures / output.name)
        print(f"Rendered {len(outputs)} SVGs into {figures}")


if __name__ == "__main__":
    main()
