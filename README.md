# MHWTracking.jl

A Julia package port of the local MATLAB marine heatwave tracking code. It provides daily spatial patches linked by temporal overlap, independent direct three-dimensional connectivity, regular four-dimensional connectivity, ECCO LLC four-dimensional connectivity, spatial smoothing, and event normalization.

## Installation and testing

The local package root is `/Volumes/new_drive/lagrangian/julia`. It contains `Project.toml`, `src/`, `test/`, `examples/`, `validation/`, and `notebooks/`, alongside the original MATLAB files and input data retained during relocation. From a terminal, enter the package directory:

```sh
cd /Volumes/new_drive/lagrangian/julia
```

Run from the package root:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.test()
using MHWTracking
```

In another project, use `Pkg.develop(path="/Volumes/new_drive/lagrangian/julia")`, or substitute the package's location on your machine. Both the package and module are named `MHWTracking`; the minimum Julia version is 1.10. MAT.jl reads v5/v7.3 MAT files. MiniQhull.jl builds the Delaunay meshes in `spn`; DelaunayTriangulation.jl supplies adaptive triangle-containment predicates.

## Real three-dimensional data

```julia
using MHWTracking
d = load_test_data(pwd();
                   binary_rule=:positive)
tracks = track_mhw_3d(d.mask, d.lon, d.lat, 4)
println(length(tracks))
println(tracks[1].day)
```

The inspected input is `mhw_ts::Array{Float64,3}` in `mhw_ts.mat`, with dimensions 400×160×731 and coordinates in `lonlatoi.mat`. It contains NaNs and continuous values rather than a binary mask. Under the user-confirmed rule, `:positive` explicitly applies `isfinite(x) && x > 0`. `:raw` preserves all input values, including NaNs that MATLAB treats as nonzero; its output must not be interpreted as properly binarized marine heatwaves. Loading does not permute dimensions or apply smoothing.

Complete example, from the package root: `julia --project=. examples/example_3d.jl`. The example and full-input validation default to the package root for `mhw_ts.mat` and `lonlatoi.mat`; `MHW_DATA_DIR` or an explicit directory argument selects another data location. Notebook 02 uses the same environment variable. Original root-level MATLAB files and `.mat` data are local inputs excluded from Git; the preserved sources in `validation/matlab` and small fixtures in `test/fixtures` remain part of the package. The argument `4` is the minimum daily patch size used for this validation and can be adjusted for the scientific application.

## Four methods and their outputs

| Function | Input order | Algorithm | Output |
|---|---|---|---|
| `track_mhw_3d(mask,lon,lat,nums=1)` | lon×lat×time | Daily 8-connectivity with periodic longitude; overlap ≥0.5; source split/merge logic | `Vector{Track3D}` |
| `track_mhw_3d_connected(mask)` | lon×lat×time | Independent direct 26-connectivity, optionally 6; periodic longitude | `(tracks,info)`, `idx3` |
| `track_mhw_4d_lonlat(mask)` | lon×lat×depth×time | 80-connectivity; periodic longitude | `(tracks,info)`, `idx4` |
| `track_mhw_4d_llc(mask,topology)` | i×j×tile×depth×time | Per-tile 80-connectivity plus exact seams and depth/time offsets ±1 | `(tracks,info)`, `idx5` |

All indices are one-based and column-major, matching MATLAB `sub2ind`. Latitude, depth, and time are not periodic. Nonzero numeric values, including NaN, are active; apply any land or validity mask explicitly when constructing the binary input. Julia arrays must retain the dimensions specified by the API, including singleton depth and time axes.

`Track3D` retains `day`, `xloc`, `yloc`, `ori_day`, `ori_order`, `split_num`, and `split_day`. Coordinates are returned as Float64 values representing the same integer coordinates as MATLAB. Events follow MATLAB's closure order. The default minimum duration is 5; `min_duration` and `overlap_ratio` can explicitly override the original fixed thresholds. MATLAB requires `nums`, whereas Julia adds a convenience default of 1. During merging, existing events remain separate and receive portions of the current patch; splitting is recorded within the existing event. This event definition differs from direct connectivity.

`IndexedTrack` retains `day::Vector{Int32}`, `idxN::Vector{Vector{UInt64}}`, and `nvoxels::Vector{UInt32}`. The `info` named tuple contains the original MATLAB statistics, including component counts, seam counts, and active boundary cells. Keywords:

- `min_duration=5`: filter by the number of unique output days.
- `time=[]`: local-to-global day labels, which may be nonconsecutive or repeated; connectivity always follows adjacent **local** time positions.
- `mask_size=size(mask)`: global dimensions used for global linear indices.
- `save_file=""`: optional MAT output containing `tracks` and `info`, using MAT.jl's file format.
- `verbose=false`: enable progress summaries.
- Regular four-dimensional tracking also supports `add_idx5_alias=false`.

```julia
mask = falses(8, 6, 3, 7)
mask[1, 3, 1, 1] = true
mask[8, 4, 2, 2] = true
tracks, info = track_mhw_4d_lonlat(mask; min_duration=1)
@assert length(tracks) == 1
```

## Full Tasman Sea normalization example

Notebook [07_tasman_sea_spn](notebooks/07_tasman_sea_spn.ipynb) reproduces the original 230-day 2015–2016 Tasman Sea case with SST anomaly, both normalization modes, all five normalized times, three full-event GIFs, and the separate original global smoothing example. The five supplied MAT files are retained under `data/tasman_sea` and excluded from Git. See the [reproduction guide](docs/tasman_sea_example.md) for data alignment, MATLAB comparisons, and cache controls.

`load_matlab_tracks(path)` reads the original MATLAB struct/cell arrays into `Vector{Track3D}` while retaining all seven fields and one-based indices. Run `julia --project=. examples/example_tasman_spn.jl` from the package root, or open notebook 07. Set `MHW_TASMAN_DATA_DIR` to use another directory containing the five files.

## LLC topology

```julia
topology = load_seam_edges("test/fixtures/llc90_topology.mat")
a, b = topology.seam_edges[1, :]
n2 = prod(topology.grid_size)
mask = SparseMask((90,90,13,2,3), [a, b+n2+n2*2])
tracks, info = track_mhw_4d_llc(mask, topology; min_duration=1)
@assert length(tracks) == 1
```

The package includes the actual LLC90 grid and its 2,312 seam pairs, with XC/YC/Z from the user-provided files. `seamEdges` is an N×2 table of horizontal `(i,j,tile)` linear indices; the explicit pairs encode orientation and endpoint mappings. Consecutive tile numbers do not imply spatial adjacency, and no additional lateral seam diagonals are introduced. Reuse `LLCTopology` across calls. `SparseMask` supports sparse large inputs and synthetic tests; ordinary arrays are also accepted.

## Auxiliary functions

```julia
indices = grid_nearby_index([0.,1.,90.,359.], [-2.,0.,3.], 2.)
smoothed = convzz(ones(4,3), indices, trues(4,3))
```

`grid_nearby_index` selects neighbors within coordinate windows ±res, preserving the original meshgrid index order and 0..360 longitude logic. `convzz` applies a neighborhood nanmean at centers whose ocean mask is nonzero; it does not additionally mask neighboring land cells.

`spn(tracks,lon,lat,anomaly,norm_flag)` returns the original five MATLAB outputs: maximum event radii, daily spatial samples, a 100×100×5×events×1 normalized array, and unit x/y grids. `anomaly` must start at the earliest day across all tracks; a 4-D input uses only its first variable. Arithmetic coordinate centers, planar distance approximation, 100×100 polar sampling, and time extrapolation from `(1:n)/n` follow the local source. There is no added longitude correction or area weighting. With `norm_flag==0`, values outside the event boundary are zero; other flags retain the field. Empty events, fewer than two time samples, and degenerate interpolation points raise explicit errors.

## Validation and performance

Seven [Julia Jupyter notebooks](notebooks/README.md) include recorded execution outputs for event evolution, differences between the three-dimensional methods, four-dimensional boundaries, LLC seams, `spn` comparisons, and performance. Browse the [HTML gallery and notebook navigation](notebooks/gallery.html). All figures are exported as PNG/PDF; plotting and IJulia use a separate environment without adding core package dependencies.

See the [source mapping](docs/algorithm_mapping.md), [validation report](docs/validation_report.md), and [performance report](docs/performance.md). Reproducible scripts:

```sh
julia --project=. test/runtests.jl
julia --project=. validation/validate_seams.jl
julia --project=. validation/validate_real.jl
julia --project=. benchmark/benchmark_tracking.jl
```

`validation/matlab/generate_reference.m` calls the preserved original source to generate small fixtures; `real_reference.m` generates the full real-data reference. Synthetic four-dimensional validation does not establish equivalence on real production four-dimensional inputs. The independent three-dimensional API is an explicitly documented 6/26-connectivity extension: the local and GitHub source contain only an embedded two-dimensional helper, not the standalone three-dimensional file described in the request.

Methods: [Sun et al. 2023](https://doi.org/10.1016/j.pocean.2022.102947), [Zhao et al. 2026](https://doi.org/10.1038/s41598-026-40354-4). Reference code: [MHW_tracking](https://github.com/ZijieZhaoMMHW/MHW_tracking). The downloaded `mhwtrack`/`spn` files match the local files. The repository README gives a different argument order from the source; this package uses the source's mask-first order.

## License

GPL-3.0, following the reference MATLAB repository; see LICENSE. Attribution to original author Sun Di and modifying author Zijie Zhao is retained.
