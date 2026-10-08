"""Execute genuine IJulia kernels, preserve outputs, and export standalone HTML."""
import argparse
import json
import os
from pathlib import Path
import shutil
import time
import runpy

import nbformat
from nbclient import NotebookClient
from nbconvert import HTMLExporter
from jupyter_client import KernelManager
from jupyter_client.kernelspec import KernelSpecManager

DIRECTORY = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("notebooks", nargs="*", help="Notebook basenames; default: all available notebooks")
    parser.add_argument("--julia", default=shutil.which("julia"), help="Julia executable")
    args = parser.parse_args()
    if not args.julia:
        parser.error("Julia executable not found; pass --julia")

    # Keep registration/runtime files inside this workspace, not in user-wide Jupyter config.
    data = DIRECTORY / ".jupyter"
    kernel_dir = data / "kernels" / "mhwtracking-local"
    kernel_dir.mkdir(parents=True, exist_ok=True)
    runtime = data / "runtime"
    runtime.mkdir(parents=True, exist_ok=True)
    os.environ["JUPYTER_RUNTIME_DIR"] = str(runtime)
    os.environ["IPYTHONDIR"] = str(data / "ipython")
    spec = {
        "argv": [args.julia, "-i", "--startup-file=no", f"--project={DIRECTORY}",
                 "-e", "import IJulia; IJulia.run_kernel()", "{connection_file}"],
        "display_name": "MHWTracking Julia (local)", "language": "julia", "env": {},
    }
    (kernel_dir / "kernel.json").write_text(json.dumps(spec, indent=2))
    paths = [DIRECTORY / name for name in args.notebooks] if args.notebooks else sorted(DIRECTORY.glob("[0-9][0-9]_*.ipynb"))
    html_dir = DIRECTORY / "html"
    html_dir.mkdir(exist_ok=True)
    report = DIRECTORY / "execution_report.json"
    previous = json.loads(report.read_text()) if report.exists() else []
    selected_names = {p.name for p in paths}
    summary = [entry for entry in previous if entry["notebook"] not in selected_names]
    for path in paths:
        print(f"Executing {path.name} ...", flush=True)
        notebook = nbformat.read(path, as_version=4)
        manager = KernelManager(kernel_name="mhwtracking-local",
            kernel_spec_manager=KernelSpecManager(kernel_dirs=[str(data / "kernels")]))
        client = NotebookClient(notebook, km=manager, timeout=900, startup_timeout=300,
            resources={"metadata": {"path": str(DIRECTORY)}}, allow_errors=False, record_timing=True)
        start = time.monotonic()
        try:
            client.execute()
        finally:
            # Preserve diagnostic outputs even on failure; failure still stops the run.
            nbformat.write(notebook, path)
            # A supplied manager is owned by us, including when a cell fails.
            if manager.has_kernel:
                manager.shutdown_kernel(now=True)
        nbformat.validate(notebook)
        errors = [o for c in notebook.cells if c.cell_type == "code" for o in c.outputs if o.output_type == "error"]
        if errors:
            raise RuntimeError(f"{path.name}: notebook errors found")
        images = sum("image/png" in o.get("data", {}) for c in notebook.cells if c.cell_type == "code" for o in c.outputs)
        html, _ = HTMLExporter(template_name="lab").from_notebook_node(notebook)
        (html_dir / f"{path.stem}.html").write_text(html)
        entry = {"notebook": path.name, "code_cells": sum(c.cell_type == "code" for c in notebook.cells),
                 "png_outputs": images, "errors": len(errors), "seconds": round(time.monotonic() - start, 2)}
        summary.append(entry)
        print(json.dumps(entry), flush=True)
    report.write_text(json.dumps(sorted(summary, key=lambda entry: entry["notebook"]), indent=2) + "\n")
    runpy.run_path(str(DIRECTORY / "build_gallery.py"))


if __name__ == "__main__":
    main()
