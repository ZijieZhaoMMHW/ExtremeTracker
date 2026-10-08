function time_options(mask,time,mask_size,min_duration)
    dims=size(mask); full=Tuple(Int.(mask_size)); n=ndims(mask)
    length(full)==n && full[1:n-1]==dims[1:n-1] || throw(DimensionMismatch("mask_size spatial dimensions must match"))
    min_duration>=1 || throw(ArgumentError("min_duration must be >= 1"))
    days=isempty(time) ? collect(1:dims[end]) : vec(collect(time))
    length(days)==dims[end] || throw(DimensionMismatch("time length must match local time"))
    all(d->isinteger(d) && 1<=d<=full[end],days) || throw(ArgumentError("time labels must be integer indices within mask_size"))
    Int.(days),full
end
function assemble_tracks(::Val{N},comps,dsu,days,nspace,min_duration;alias=false) where N
    roots=[root!(dsu,c) for c in eachindex(comps)]
    rootids=sort!(unique(roots)); lookup=Dict(r=>i for (i,r) in enumerate(rootids))
    buckets=[Dict{Int,Vector{UInt64}}() for _ in rootids]
    for (c,pixels) in enumerate(comps)
        bucket=buckets[lookup[roots[c]]]
        for p in pixels
            localday=fld(p-1,nspace)+1; gd=days[localday]
            idx=UInt64(mod1(p,nspace))+UInt64(gd-1)*UInt64(nspace)
            push!(get!(bucket,gd,UInt64[]),idx)
        end
    end
    tracks=IndexedTrack{N}[]
    for bucket in buckets
        d=sort!(collect(keys(bucket))); length(d)>=min_duration || continue
        ix=[bucket[day] for day in d]
        push!(tracks,IndexedTrack{N}(Int32.(d),ix,UInt32.(length.(ix)),alias))
    end
    tracks,length(rootids)
end
function component_map(cc)
    map=zeros(Int32,cc.ImageSize)
    for (c,p) in enumerate(cc.PixelIdxList); map[p].=c; end
    map
end
function apply_pairs!(dsu,pairs,n)
    ps=sort!(collect(pairs);by=p->p[1]+n*(p[2]-1))
    for (a,b) in ps; union!(dsu,a,b); end
    length(ps)
end

"""
    track_mhw_4d_lonlat(mask; time=[], mask_size=size(mask), min_duration=5,
                        add_idx5_alias=false, save_file="", verbose=false)

Full 80-connected `[lon,lat,depth,time]` components with periodic longitude,
bounded latitude/depth/time. Time labels affect global idx4 and occupied-day
filtering, not local adjacency. Returns `(tracks, info)`; exact MATLAB root
ordering and per-day concatenation are retained. NaN is nonzero and active.
"""
function track_mhw_4d_lonlat(mask::AbstractArray{<:Number,4};time=Int[],
        mask_size=size(mask),min_duration::Real=5,add_idx5_alias=false,save_file="",verbose=false)
    nx,ny,nz,nt=size(mask); nx>=2 || throw(ArgumentError("longitude needs at least two cells"))
    isempty(mask) && throw(ArgumentError("mask must be nonempty"))
    days,full=time_options(mask,time,mask_size,min_duration)
    cc=components(mask); map=component_map(cc); d=DSU(cc.NumObjects)
    pairs=Set{Tuple{Int,Int}}(); active_boundary=0
    for t in 1:nt,z in 1:nz,y in 1:ny,x in (1,nx)
        a=Int(map[x,y,z,t]); a==0 && continue; active_boundary+=1
        xx=x==1 ? nx : 1
        for dt in -1:1,dz in -1:1,dy in -1:1
            yy=y+dy; zz=z+dz; tt=t+dt
            (1<=yy<=ny && 1<=zz<=nz && 1<=tt<=nt) || continue
            b=Int(map[xx,yy,zz,tt]); b>0 && b!=a && push!(pairs,(a,b))
        end
    end
    npairs=apply_pairs!(d,pairs,cc.NumObjects)
    tracks,nraw=assemble_tracks(Val(4),cc.PixelIdxList,d,days,nx*ny*nz,min_duration;alias=Bool(add_idx5_alias))
    info=(method="regular_lonlat_4d_bwconncomp_periodic_lon_union",mask_size=full,
        grid_size=(nx,ny,nz),ntime_local=nt,ntime_total=full[end],min_duration=min_duration,
        connectivity="full 3x3x3x3 within regular lon-lat grid; periodic lon seam union with dlat,dz,dt in -1:1",
        index_field="idx4",has_idx5_alias=Bool(add_idx5_alias),
        n_components_before_periodic_union=cc.NumObjects,n_raw_components=nraw,n_tracks=length(tracks),
        n_lon_seam_edges=2sum(length(max(1,y-1):min(ny,y+1)) for y in 1:ny),
        n_active_lon_boundary_cells=active_boundary,n_lon_seam_component_pairs=npairs)
    verbose && println("Regular 4-D: $(cc.NumObjects) local components, $nraw unioned, $(length(tracks)) tracks")
    isempty(save_file) || save_tracks(save_file,tracks,info)
    tracks,info
end

"""
    track_mhw_4d_llc(mask, topology_or_seam_edges; time=[], mask_size=size(mask),
                     min_duration=5, save_file="", verbose=false)

Input `[i,j,tile,depth,time]`. Full 80 connectivity within each tile, then exact
seam edges with dz,dt ∈ -1:1. No ordinary adjacency along the tile axis. Reuse
an LLCTopology to avoid rebuilding and validating the horizontal topology.
Returns `(tracks, info)` with one-based global `idx5` per occupied global day.
"""
function track_mhw_4d_llc(mask::AbstractArray{<:Number,5},topology;
        time=Int[],mask_size=size(mask),min_duration::Real=5,save_file="",verbose=false)
    isempty(mask) && throw(ArgumentError("mask must be nonempty"))
    ni,nj,tiles,nz,nt=size(mask); n2=ni*nj*tiles; nspace=n2*nz
    topo=topology isa LLCTopology ? topology : LLCTopology(topology,(ni,nj,tiles))
    topo.grid_size==(ni,nj,tiles) || throw(DimensionMismatch("topology does not match mask"))
    days,full=time_options(mask,time,mask_size,min_duration)
    comps=Vector{Vector{Int}}(); boundary=Dict{Int,Int}(); byhorizontal=Dict{Int,Vector{Int}}()
    tilepixels=[Int[] for _ in 1:tiles]; ci5=CartesianIndices(size(mask))
    tile_dims=(ni,nj,nz,nt); tile_li=LinearIndices(tile_dims)
    for p in active_indices(mask)
        i,j,tile,z,t=Tuple(ci5[p]); push!(tilepixels[tile],tile_li[i,j,z,t])
    end
    for tile in 1:tiles
        cc=components(SparseMask(tile_dims,tilepixels[tile])); ci=CartesianIndices(cc.ImageSize)
        for pixels in cc.PixelIdxList
            id=length(comps)+1; globalpixels=Int[]
            for p in pixels
                i,j,z,t=Tuple(ci[p]); h=i+(j-1)*ni+(tile-1)*ni*nj
                gp=h+(z-1)*n2+(t-1)*nspace; push!(globalpixels,gp)
                if haskey(topo.neighbors,h)
                    boundary[gp]=id; push!(get!(byhorizontal,h,Int[]),gp)
                end
            end
            push!(comps,globalpixels)
        end
    end
    d=DSU(length(comps)); npairs=0
    for row in eachrow(topo.seam_edges)
        a,b=row; pairs=Set{Tuple{Int,Int}}()
        for p in get(byhorizontal,a,Int[])
            z=mod(fld(p-1,n2),nz)+1; t=fld(p-1,nspace)+1; source=boundary[p]
            for dt in -1:1,dz in -1:1
                zz=z+dz; tt=t+dt; 1<=zz<=nz && 1<=tt<=nt || continue
                target=get(boundary,b+(zz-1)*n2+(tt-1)*nspace,0)
                target>0 && target!=source && push!(pairs,(source,target))
            end
        end
        npairs+=apply_pairs!(d,pairs,length(comps))
    end
    tracks,nraw=assemble_tracks(Val(5),comps,d,days,nspace,min_duration)
    info=(method="tilewise_4d_bwconncomp_exact_llc_seam_union",mask_size=full,
        grid_size=(ni,nj,tiles,nz),ntime_local=nt,ntime_total=full[end],min_duration=min_duration,
        connectivity="full 3x3x3x3 within tile; exact seam union with dz,dt in -1:1",
        n_tile_components=length(comps),n_raw_components=nraw,n_tracks=length(tracks),
        n_seam_edges=size(topo.seam_edges,1),n_active_seam_cells=length(boundary),n_seam_component_pairs=npairs)
    verbose && println("LLC 4-D: $(length(comps)) tile components, $nraw unioned, $(length(tracks)) tracks")
    isempty(save_file) || save_tracks(save_file,tracks,info)
    tracks,info
end

"""
    track_mhw_3d_connected(mask; connectivity=26, min_duration=5, time=[], ...)

Independent direct space-time method, with 6 or 26 connectivity and periodic
longitude. It intentionally differs from daily overlap tracking. This is a
documented extension because no independent MATLAB 3-D helper was supplied.
Returns `(tracks, info)` exposing per-day global idx3.
"""
function track_mhw_3d_connected(mask::AbstractArray{<:Number,3};connectivity=26,
        min_duration::Real=5,time=Int[],mask_size=size(mask),save_file="",verbose=false)
    days,full=time_options(mask,time,mask_size,min_duration)
    cc=bwconncomp_periodic_lon(mask;connectivity)
    tracks,nraw=assemble_tracks(Val(3),cc.PixelIdxList,DSU(cc.NumObjects),days,prod(size(mask)[1:2]),min_duration)
    info=(method="direct_3d_periodic_lon_extension",mask_size=full,connectivity=connectivity,
        n_raw_components=nraw,n_tracks=length(tracks),min_duration=min_duration)
    verbose && println("Direct 3-D: $nraw components, $(length(tracks)) tracks")
    isempty(save_file) || save_tracks(save_file,tracks,info)
    tracks,info
end
