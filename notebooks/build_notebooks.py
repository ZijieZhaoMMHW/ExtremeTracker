"""Build Julia notebooks from reviewable paired .jl cell sources."""
from pathlib import Path
import hashlib
import nbformat

DIRECTORY = Path(__file__).resolve().parent


def build(path):
    cells, lines, kind = [], [], "code"

    def flush():
        source = "\n".join(lines).strip()
        if not source:
            return
        factory = nbformat.v4.new_markdown_cell if kind == "markdown" else nbformat.v4.new_code_cell
        cell = factory(source)
        cell.id = hashlib.sha256(f"{path.stem}:{len(cells)}:{source}".encode()).hexdigest()[:16]
        cells.append(cell)

    for line in path.read_text().splitlines():
        if line.startswith("# %%"):
            flush()
            lines = []
            kind = "markdown" if "[markdown]" in line else "code"
        elif kind == "markdown":
            lines.append(line[2:] if line.startswith("# ") else line.removeprefix("#"))
        else:
            lines.append(line)
    flush()
    notebook = nbformat.v4.new_notebook(cells=cells, metadata={
        "kernelspec": {"display_name": "Julia 1.12", "language": "julia", "name": "julia-1.12"},
        "language_info": {"name": "julia", "version": "1.12.3", "file_extension": ".jl", "mimetype": "application/julia"},
        "mhwtracking": {"source": path.name, "documentation_language": "en",
                        "description": "Generated from the paired Julia cell source; executed with IJulia."},
    })
    output = path.with_suffix(".ipynb")
    nbformat.validate(notebook)
    nbformat.write(notebook, output)
    print(f"Built {output.name}: {len(cells)} cells", flush=True)


if __name__ == "__main__":
    for path in sorted(DIRECTORY.glob("[0-9][0-9]_*.jl")):
        build(path)
