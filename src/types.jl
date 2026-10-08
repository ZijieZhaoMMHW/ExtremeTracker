"""MATLAB-compatible overlap track, with one coordinate vector per occupied day."""
struct Track3D
    day::Vector{Int}
    xloc::Vector{Vector{Float64}}
    yloc::Vector{Vector{Float64}}
    ori_day::Int
    ori_order::Int
    split_num::Vector{Int}
    split_day::Vector{Int}
end

"""Connected event. `idx3`, `idx4`, or `idx5` exposes one-based global indices per day."""
struct IndexedTrack{N}
    day::Vector{Int32}
    indices::Vector{Vector{UInt64}}
    nvoxels::Vector{UInt32}
    has_idx5_alias::Bool
end
function Base.getproperty(t::IndexedTrack{N}, s::Symbol) where N
    if s == Symbol("idx", N) || (s == :idx5 && getfield(t, :has_idx5_alias))
        return getfield(t, :indices)
    end
    getfield(t, s)
end
Base.propertynames(t::IndexedTrack{N}) where N =
    t.has_idx5_alias ? (:day, Symbol("idx", N), :nvoxels, :idx5) : (:day, Symbol("idx", N), :nvoxels)

"""Connected components in MATLAB column-major order; indices are one-based."""
struct ConnectedComponents{N}
    Connectivity::Int
    ImageSize::NTuple{N,Int}
    NumObjects::Int
    PixelIdxList::Vector{Vector{Int}}
end

mutable struct DSU
    parent::Vector{Int}
    rank::Vector{UInt8}
end
DSU(n) = DSU(collect(1:n), zeros(UInt8,n))
function root!(d::DSU, x::Int)
    r=x
    while d.parent[r] != r
        r=d.parent[r]
    end
    while d.parent[x] != x
        nxt=d.parent[x]; d.parent[x]=r; x=nxt
    end
    r
end
function union!(d::DSU,a::Int,b::Int)
    a=root!(d,a); b=root!(d,b)
    a==b && return
    if d.rank[a]<d.rank[b]
        d.parent[a]=b
    elseif d.rank[a]>d.rank[b]
        d.parent[b]=a
    else
        d.parent[b]=a; d.rank[a]+=1
    end
end

"""Canonical full voxel membership, independent of event labels."""
function event_partition(tracks::AbstractVector{Track3D}, dims)
    nx,ny=dims[1:2]; n=nx*ny
    parts=[sort!([Int(t.xloc[j][k])+(Int(t.yloc[j][k])-1)*nx+(t.day[j]-1)*n
          for j in eachindex(t.day) for k in eachindex(t.xloc[j])]) for t in tracks]
    sort!(parts)
end
event_partition(tracks::AbstractVector{<:IndexedTrack}, dims=nothing) =
    sort!([sort!(reduce(vcat,t.indices;init=UInt64[])) for t in tracks])
"""Compare exact event membership, not just component counts."""
partition_equal(a,b,dims=nothing) = event_partition(a,dims)==event_partition(b,dims)

"""Build a dense label array. Throws if two tracks claim the same voxel."""
function labels(tracks, dims)
    out=zeros(Int32,Tuple(dims))
    for (k,t) in enumerate(tracks)
        indices=t isa Track3D ? event_partition([t],dims)[1] : reduce(vcat,t.indices;init=UInt64[])
        for p in indices
            out[p]==0 || throw(ArgumentError("overlapping event memberships at index $p"))
            out[p]=k
        end
    end
    out
end
