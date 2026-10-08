# Julia Jupyter notebooks and visualizations

The seven `.ipynb` files contain recorded IJulia execution outputs, numeric validation assertions, and embedded figures. Paired `.jl` files provide reviewable cell sources; the original example and validation scripts remain independently runnable. Every figure is also exported as PNG and PDF for reports and manuscript drafts.

Execution on 2026-10-08: **seven notebooks, all 40 code cells executed successfully, zero errors; 21 embedded PNG outputs, 21 PNG/PDF figure pairs, seven standalone HTML exports, and three 230-frame GIFs**. The real three-dimensional notebook again matched all 4,277 MATLAB events exactly; tracking took approximately 3.10 seconds in that run. Versions are recorded in `environment_versions.json`. Notebook execution time includes loading, compilation, and plotting and is not the tracker-only time.

| Notebook | Contents | External input |
|---|---|---|
| [01_tracking_concepts.ipynb](01_tracking_concepts.ipynb) | Periodic boundaries, split/merge behavior, differences between the two 3-D methods | Packaged MATLAB fixtures |
| [02_real_3d.ipynb](02_real_3d.ipynb) | Full real input, occurrence frequency, event statistics/evolution, full MATLAB comparison | `mhw_ts.mat` / `lonlatoi.mat`; full MATLAB output is optional |
| [03_regular_4d.ipynb](03_regular_4d.ipynb) | 80-connectivity, longitude/depth/time boundaries, global-time mapping | Packaged MATLAB fixtures |
| [04_llc_seams.ipynb](04_llc_seams.ipynb) | LLC geography, tile-pair counts, orientation, all 11,560 seam cases | Packaged actual grid/seams and MATLAB fixtures |
| [05_auxiliary_spn.ipynb](05_auxiliary_spn.ipynb) | Neighborhood smoothing, real events under both norm_flag values, five normalized times, MATLAB residuals | Packaged MATLAB fixtures |
| [06_validation_performance.ipynb](06_validation_performance.ipynb) | Recomputed strict comparisons, recorded full-input timing, current synthetic 4-D measurements | Large recorded metrics are optional; small tests need no external input |
| [07_tasman_sea_spn.ipynb](07_tasman_sea_spn.ipynb) | Full 230-day Tasman Sea event, both spatial modes, normalized times, full MATLAB checks, five figures and three GIFs | Five original MAT files under `data/tasman_sea`; full MATLAB references are optional |

## Complete original Tasman Sea case

[07_tasman_sea_spn.ipynb](07_tasman_sea_spn.ipynb) adds the full 2015-10-03 to 2016-05-19 event: 230 days of supplied SST anomaly, both spatial modes, five normalized times, geographic evolution, three complete GIFs, and the separate global smoothing example. It uses all five supplied MAT files in `data/tasman_sea`; set `MHW_TASMAN_DATA_DIR` for another location. Full outputs use an input/code SHA256 cache; `MHW_TASMAN_RECOMPUTE=1` forces recomputation. See [the reproduction guide](../docs/tasman_sea_example.md). The original six notebooks retain their recorded execution outputs. Notebook 07 executed all nine code cells without errors; both complete normalization modes matched the original MATLAB outputs, and the separate smoothing checks passed. `animation_qa.json` records 230 frames at 500 ms/frame for every GIF.

## Browse the results

- [HTML gallery and notebook navigation](gallery.html): links to seven notebooks, their individual figures, and the full Tasman Sea animations.
- Each file in `html/` embeds its figures and can be shared independently with readers who do not have Julia/Jupyter.
- `figures/` contains standalone PNG and PDF exports. `execution_report.json` records code-cell counts, PNG outputs, and errors from actual execution.

## Run interactively

From the package root `/Volumes/new_drive/lagrangian/julia` (Julia ≥1.10; locally tested with 1.12.3):

```sh
cd /Volumes/new_drive/lagrangian/julia
# Use the installed workspace depot while retaining the user's existing depot
export JULIA_DEPOT_PATH="$PWD/.julia_depot:$HOME/.julia"
julia notebooks/setup.jl
python3 -m pip install -r notebooks/requirements.txt  # Only if Jupyter is not installed
jupyter lab notebooks
```

Select an installed **Julia kernel** in Jupyter, open a notebook, and choose Run All. The first cell activates the separate `notebooks/Project.toml` environment. If no Julia kernel is installed, run in a terminal:

```sh
julia --project=notebooks -e 'using IJulia; IJulia.installkernel("MHWTracking Julia", "--project=@.")'
```

Notebook 02 defaults to `mhw_ts.mat` and `lonlatoi.mat` in the package root, alongside the original MATLAB files retained during relocation. Change `DATA_DIR` in its configuration cell or set `MHW_DATA_DIR` before starting Jupyter to use another data directory. Saved outputs can be viewed without the external data; notebooks 01 and 03–06 can be rerun using packaged fixtures; notebook 07 requires the five Tasman Sea input files. Cells requiring absent full MATLAB outputs or original timing files explicitly print NOT RUN.

The plotting library incurs compilation overhead on first load. Notebook 02 reloads the approximately 374 MB raw array and retracks the full input; loading and plotting are excluded from tracker timing. Array/event statistics are not area-weighted, and day indices are not converted to unknown calendar dates.

## Rebuild, execute, and export

```sh
python3 notebooks/build_notebooks.py
python3 notebooks/execute_notebooks.py --julia /absolute/path/to/julia
# Debug a single notebook, for example:
python3 notebooks/execute_notebooks.py 05_auxiliary_spn.ipynb --julia /absolute/path/to/julia
```

Rebuilding clears previous outputs; execute afterward to update results. The batch runner uses a workspace-local temporary kernelspec without changing global Jupyter configuration. It starts actual IJulia kernels and stops on errors; outputs are not simulated. `execution_report.json` retains the latest selected notebook records and preserves records for unselected notebooks. Execution updates the gallery automatically; rebuild only the gallery with `python3 notebooks/build_gallery.py`.

**CairoMakie** and **IJulia** belong to the separate notebook environment and are not core package dependencies. See the [official IJulia installation guide](https://ijulia.org/stable/manual/installation/) and [CairoMakie documentation](https://docs.makie.org/stable/explanations/backends/cairomakie).

## Scientific and validation scope

The real mask uses the confirmed rule `isfinite(x) && x > 0`. Real 3-D tracking preserves MATLAB's first/last longitude connection even though the input region is not global; figures and explanations retain this fact. Direct 3-D connectivity is a documented extension. Four-dimensional MHW figures use synthetic cases; actual LLC grid/seam pairs come from the user's data. Notebook 05's `spn` numeric field is the original mhw_ts and is not assigned an unconfirmed SST-anomaly interpretation. Notebook 07 uses the author's actual Tasman Sea SST anomaly in °C, with the original calendar-date convention. Recorded timings and current measurements are identified separately; cumulative allocation and maximum RSS measure different quantities.

More details: [validation report](../docs/validation_report.md), [performance report](../docs/performance.md).
