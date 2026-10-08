# Nonzero includes NaN, matching mask ~= 0 / bwconncomp numeric input.
active(x::Number) = x != 0
active(x) = throw(ArgumentError("masks must contain numeric or Boolean values"))

"""
    SparseMask(dims, active_linear_indices)

Read-only binary array storing only active one-based linear indices. Useful for
large, sparse LLC inputs and isolated seam verification; preserves array layout.
"""
struct SparseMask{N} <: AbstractArray{Bool,N}
    dims::NTuple{N,Int}
    indices::Vector{Int}
    lookup::Set{Int}
end
function SparseMask(dims::NTuple{N,Int},indices) where N
    all(>=(0),dims) || throw(ArgumentError("nonnegative dimensions required"))
    ps=sort!(unique(Int.(collect(indices))))
    all(p->1<=p<=prod(dims),ps) || throw(BoundsError(dims,ps))
    SparseMask(dims,ps,Set(ps))
end
Base.size(a::SparseMask)=a.dims
Base.IndexStyle(::Type{<:SparseMask})=IndexLinear()
Base.getindex(a::SparseMask,p::Int)=(checkbounds(a,p); p in a.lookup)
active_indices(a::SparseMask)=a.indices
active_indices(a::BitArray)=findall(vec(a))
active_indices(a::AbstractArray)=findall(active,vec(a))

@inline function within_bounds(q,dims)
    for d in eachindex(dims)
        1<=q[d]<=dims[d] || return false
    end
    true
end

function components(mask::SparseMask{N};periodic=false,connectivity=3^N-1) where N
    connectivity in (2N,3^N-1) || throw(ArgumentError("unsupported connectivity"))
    offsets=[Tuple(c) for c in CartesianIndices(ntuple(_->-1:1,N))
        if any(!iszero,Tuple(c)) && (connectivity==3^N-1 || sum(abs,Tuple(c))==1)]
    ci=CartesianIndices(size(mask)); li=LinearIndices(size(mask)); seen=Set{Int}(); parts=Vector{Vector{Int}}()
    for seed in mask.indices
        seed in seen && continue
        push!(seen,seed); queue=Int[seed]; head=1
        while head<=length(queue)
            p=queue[head]; head+=1; coord=Tuple(ci[p])
            for off in offsets
                q=ntuple(d->coord[d]+off[d],N)
                periodic && (q=(mod1(q[1],size(mask,1)),q[2:end]...))
                within_bounds(q,size(mask)) || continue
                ix=li[q...]
                if ix ∉ seen && ix in mask.lookup
                    push!(seen,ix); push!(queue,ix)
                end
            end
        end
        push!(parts,sort!(queue))
    end
    ConnectedComponents(connectivity,size(mask),length(parts),parts)
end

function components(mask::AbstractArray{T,N}; periodic=false, connectivity=3^N-1) where {T,N}
    connectivity in (2N,3^N-1) || throw(ArgumentError("use axial $(2N) or full $(3^N-1) connectivity"))
    dims=size(mask); ci=CartesianIndices(mask); li=LinearIndices(mask)
    offsets=[Tuple(c) for c in CartesianIndices(ntuple(_->-1:1,N))
        if any(!iszero,Tuple(c)) && (connectivity==3^N-1 || sum(abs,Tuple(c))==1)]
    seen=falses(dims); parts=Vector{Vector{Int}}()
    for seed in eachindex(seen)
        (seen[seed] || !active(mask[seed])) && continue
        seen[seed]=true; queue=Int[seed]; head=1
        while head<=length(queue)
            p=queue[head]; head+=1; coord=Tuple(ci[p])
            for off in offsets
                q=ntuple(d->coord[d]+off[d],N)
                if periodic
                    q=(mod1(q[1],dims[1]),q[2:end]...)
                end
                within_bounds(q,dims) || continue
                ix=li[q...]
                if !seen[ix] && active(mask[ix])
                    seen[ix]=true; push!(queue,ix)
                end
            end
        end
        push!(parts,sort!(queue))
    end
    ConnectedComponents(connectivity,dims,length(parts),parts)
end

"""
    bwconncomp_periodic_lon(mask; connectivity=ndims(mask)==2 ? 8 : 26)

For 2-D input, reproduce the local MATLAB helper: ordinary 8-connected
components followed by periodic longitude unions, retaining concatenation order.
For 3-D input, use explicit 6/26 space-time connectivity and longitude wrapping.
The 3-D behavior is an extension: the supplied helper is only 2-D.
"""
function bwconncomp_periodic_lon(mask::AbstractArray; connectivity=ndims(mask)==2 ? 8 : 26)
    ndims(mask) in (2,3) || throw(ArgumentError("expected a 2-D or 3-D mask"))
    if ndims(mask)==3
        return components(mask;periodic=true,connectivity)
    end
    connectivity==8 || throw(ArgumentError("MATLAB daily helper uses 8 connectivity"))
    cc=components(mask); nx,ny=size(mask)
    (nx<2 || cc.NumObjects<2) && return cc
    map=zeros(Int32,size(mask))
    for (c,ps) in enumerate(cc.PixelIdxList); map[ps].=c; end
    d=DSU(cc.NumObjects)
    for y in 1:ny
        a=Int(map[1,y]); a==0 && continue
        for yy in max(1,y-1):min(ny,y+1)
            b=Int(map[nx,yy]); b==0 && continue
            ra=root!(d,a); rb=root!(d,b); d.parent[rb]=ra
        end
    end
    out=Vector{Vector{Int}}(); ids=Dict{Int,Int}()
    for (c,ps) in enumerate(cc.PixelIdxList)
        r=root!(d,c)
        if !haskey(ids,r)
            push!(out,copy(ps)); ids[r]=length(out)
        else
            append!(out[ids[r]],ps)
        end
    end
    ConnectedComponents(8,size(mask),length(out),out)
end
