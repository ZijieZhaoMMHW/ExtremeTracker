using Pkg
Pkg.activate(@__DIR__)
using MHWTracking, CairoMakie, MAT, Statistics, Test, Printf
const ROOT = normpath(joinpath(@__DIR__, ".."))
const FIGURES = joinpath(@__DIR__, "figures")
mkpath(FIGURES)
CairoMakie.activate!(type="png")
set_theme!(Theme(fontsize=15, backgroundcolor=:white,
    Axis=(xgridvisible=false, ygridvisible=false, titlesize=16,)))
include(joinpath(ROOT, "test", "reference_helpers.jl"))
include(joinpath(@__DIR__, "helpers.jl"))
println("Julia ", VERSION, " | MHWTracking ", pkgversion(MHWTracking))

