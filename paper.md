---
title: 'MHWTracking.jl: Marine heatwave tracking on regular and LLC ocean grids'
tags:
  - Julia
  - oceanography
  - marine heatwaves
  - spatiotemporal tracking
  - ocean model grids
  - reproducible research
authors:
  - name: Author name to be supplied
    affiliation: '1'
affiliations:
  - name: Affiliation and address to be supplied
    index: 1
date: 8 October 2026
bibliography: paper.bib
---

# Summary

Marine heatwaves are prolonged episodes of unusually warm ocean water [@hobday2016]. Their geographic footprints can expand, move, split, and merge, making it difficult to describe a complete event using temperature records at individual locations. MHWTracking.jl is a Julia package that groups previously identified heatwave cells into evolving ocean events and transforms surface temperature fields into a common frame for comparison. It provides surface tracking, connectivity through longitude, latitude, depth, and time, and explicit connections across the tiled latitude-longitude-cap (LLC) grids used in ocean models. Spatial smoothing, MATLAB data exchange, and event normalization connect these tracking methods to existing research workflows. Executed Jupyter notebooks demonstrate the methods, visualize their different event definitions, and reproduce a complete 230-day Tasman Sea example.

# Statement of need

Researchers studying the evolution and drivers of marine heatwaves need to follow coherent structures across space and time. A collection of independently detected temperature extremes does not itself specify which cells belong to the same moving event. Spatial overlaps, changes in shape, longitude boundaries, and subsurface connections all affect that assignment. Tiled ocean grids add another difficulty: neighboring cells can occupy different storage tiles, whose numerical ordering does not establish geographic adjacency.

MHWTracking.jl addresses these tasks for oceanographers and climate researchers working with observational fields and ocean-model output. It makes established MATLAB tracking and normalization workflows available within Julia [@bezanson2017], while retaining explicit, testable relationships to their reference implementations. The package consumes event masks; climatology estimation and temperature-threshold detection remain upstream operations. This separation lets researchers compare tracking definitions using the same input detections, or combine the package with their preferred heatwave detector. Typed results, reusable topology objects, and portable MATLAB-generated test fixtures support repeated analyses without requiring MATLAB for routine package use.

# State of the field

The Python module marineHeatWaves [@oliverCode], heatwaveR [@schlegel2018], and the MATLAB marine heatwave toolbox [@zhao2019] implement temperature-based event detection and associated statistics. These capabilities supply the local detections from which spatial events can be constructed. Sun et al. [-@sun2023] developed a framework for tracking spatially coherent marine heatwaves; Zhao et al. [-@zhao2026] introduced normalization for comparing events with different spatial scales and lifetimes. The MHW_tracking MATLAB repository provides the direct reference workflow [@mhwTrackingCode].

MHWTracking.jl complements these tools by integrating overlap tracking, direct connectivity, explicit LLC topology, and normalization in one Julia package. A separate package is justified by the need to preserve the reference workflow's event bookkeeping while supporting Julia arrays and depth-resolved tiled grids. Extending a temperature detector alone would not establish these tracking semantics. Established Julia geometry and MAT-file libraries are reused for triangulation and data exchange. The software contribution is the validated implementation, grid-aware tracking architecture, and reproducible analysis workflow; the published tracking and normalization methods retain their original attribution.

# Software design

## Event definitions and topology

The surface overlap tracker identifies daily patches using eight-neighbor connectivity with periodic longitude. Successive patches are linked when their intersection reaches the specified fraction of the smaller patch's cell count, with a default of 0.5. The original split and merge bookkeeping is retained in typed tracks, including daily membership and origin and split records. This prioritizes reproducibility of existing analyses over silently changing legacy behavior.

A separate surface method constructs connected components directly through space and time, with six- or 26-neighbor connectivity. This is an explicit extension, with a different event definition from overlap tracking. Depth-resolved tracking uses 80-neighbor connectivity in four physical dimensions. Regular grids connect their longitude endpoints; LLC inputs are processed within tiles and then united across supplied horizontal seam pairs, including adjacent depth and time positions. Separating connectivity from storage layout avoids interpreting tile numbers as geographic neighbors and supports LLC grids used by ECCO [@forget2015]. Dense masks and sparse lists of active cells share the tracking interface, while reusable seam topology avoids rebuilding grid connections for every analysis.

## Surface normalization and interoperability

Surface event normalization translates daily event fields to their arithmetic coordinate centers and scales them by the largest radius attained over the complete event. Linear scattered interpolation samples a 100-by-100 polar grid, followed by interpolation to five normalized times. One mode sets values outside the compact event boundary to zero; another retains surrounding temperature anomalies. Spatial bins accelerate triangle and boundary queries while preserving the reference interpolation and masking rules. The implementation also reproduces the source's planar distance approximation and extrapolation to normalized time zero. These conventions are documented because changing them would change composite analyses.

MATLAB import and export preserve one-based indices, array orientation, and track fields. Plotting and notebook dependencies occupy a separate environment, keeping the core package focused on tracking and normalization. Small reference fixtures support ordinary Julia package tests; larger validation drivers and notebooks support full-data reproduction. The package uses standard Julia project metadata and is distributed under GPL-3.0, following the MATLAB reference repository.

# Research impact statement

The demonstrated application is the migration and reproduction of established marine heatwave analyses. On a 400-by-160-by-731 input, converted to a mask using finite positive values, the surface overlap tracker reproduced 4,277 MATLAB events exactly, including all seven track fields, event ordering, and daily pixel ordering. The standard Julia test suite passed 716 assertions. Additional LLC validation covered 11,560 two-cell configurations across all 2,312 seams of an actual LLC90 grid. These tests provide concrete evidence that storage boundaries and implementation changes preserve the specified event assignments.

The full Tasman Sea reproduction imports the original track and sea-surface-temperature anomaly data for 3 October 2015 to 19 May 2016. Both normalization modes were compared with freshly executed original MATLAB code across every day and all five normalized times. The largest absolute finite-value difference was $1.07\times10^{-13}$ degrees Celsius, with identical NaN and zero locations. Figure 1 illustrates how the workflow converts geographic event evolution into comparable spatial and temporal coordinates.

![Reproduction of the 230-day Tasman Sea event. Upper panels show the original SST anomaly field and supplied event outline on the first, middle, and last selected days. Lower panels show event-only normalization at times 0, 0.5, and 1, using the complete event's maximum radius of 3433.24 km. Normalized time zero is extrapolated from the source time nodes. Gray denotes missing values; values outside the event boundary are zero in the lower panels. Color limits affect display only.](figures/tasman_normalization.png){width="100%"}

Measured Julia normalization calls took 27.05 and 25.86 seconds for the two modes, excluding file loading and plotting. These are single-run workflow measurements on Julia 1.12.3 running natively on ARM64; MATLAB R2020a comparisons used x86_64 compatibility execution. They are not architecture-controlled performance comparisons. Seven executed notebooks, exported figures, and MATLAB reference drivers make the demonstrated analyses inspectable. Real depth-resolved heatwave volumes and external adoption of this Julia implementation have not yet been evaluated; four-dimensional validation currently combines synthetic reference cases with the actual LLC grid and its seams.

# Acknowledgements

The package retains attribution to Di Sun and Zijie Zhao for the original and modified MATLAB tracking code. Funding acknowledgements and any additional contributor credits will be completed by the authors before submission.

# AI usage disclosure


# References
