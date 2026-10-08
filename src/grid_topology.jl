"""
    LLCTopology(seam_edges, grid_size)

Validated, reusable topology from exact one-based `(i,j,tile)` seam indices.
Orientation and corner mappings reside in the supplied pairs; no implicit tile
adjacency is added. Input row order is retained for MATLAB union ordering.
"""
struct LLCTopology
    grid_size::NTuple{3,Int}
    seam_edges::Matrix{Int}
    neighbors::Dict{Int,Vector{Int}}
end
function LLCTopology(edges::AbstractMatrix,grid_size=(90,90,13))
    dims=Tuple(Int.(grid_size)); length(dims)==3 || throw(ArgumentError("grid_size must have 3 dimensions"))
    all(>(0),dims) || throw(ArgumentError("positive grid dimensions required"))
    size(edges,2)==2 || throw(ArgumentError("seamEdges must be N×2"))
    all(x->isfinite(x) && isinteger(x) && 1<=x<=prod(dims),edges) ||
        throw(ArgumentError("seam indices must be finite integers in the horizontal grid"))
    es=Int.(edges); seen=Set{Tuple{Int,Int}}(); nbr=Dict{Int,Vector{Int}}(); ci=CartesianIndices(dims)
    for row in eachrow(es)
        a,b=row; a!=b || throw(ArgumentError("self seam edge"))
        ci[a][3]!=ci[b][3] || throw(ArgumentError("seam must cross tiles"))
        pair=minmax(a,b); pair ∉ seen || throw(ArgumentError("duplicate undirected seam edge")); push!(seen,pair)
        push!(get!(nbr,a,Int[]),b); push!(get!(nbr,b,Int[]),a)
    end
    LLCTopology(dims,es,nbr)
end

"""Read actual MAT v5/v7.3 seamEdges; return a reusable validated topology."""
load_seam_edges(path; grid_size=(90,90,13)) = LLCTopology(matread(path)["seamEdges"],grid_size)
