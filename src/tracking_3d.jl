mutable struct OpenTrack
    day::Vector{Int}
    pixels::Vector{Vector{Int}}
    single::Vector{Bool}
    ori_day::Int
    ori_order::Int
    split_num::Vector{Int}
    split_day::Vector{Int}
end
OpenTrack(day,p,order)=OpenTrack([day],[copy(p)],[true],day,order,Int[],Int[])
xy(p,nx)=(mod1.(p,nx),fld.(p.-1,nx).+1)
function output_track(t,nx)
    coords=[xy(p,nx) for p in t.pixels]
    Track3D(t.day,[Float64.(c[1]) for c in coords],[Float64.(c[2]) for c in coords],
        t.ori_day,t.ori_order,t.split_num,t.split_day)
end
function mean_lon_index_periodic(p,nx)
    theta=2pi.*(mod1.(p,nx).-1)./nx
    x=mean(cos.(theta)); y=mean(sin.(theta))
    x==0 && y==0 && return mean(mod1.(p,nx))
    a=atan(y,x); a<0 && (a+=2pi)
    idx=1+a/(2pi)*nx; idx>nx ? idx-nx : idx
end
function interp_lon_periodic(lon,idx)
    n=length(lon); f=floor(Int,idx); frac=idx-f
    l1=Float64(lon[mod1(f,n)]); l2=Float64(lon[mod1(f+1,n)])
    delta=l2-l1
    delta>180 && (l2-=360); delta< -180 && (l2+=360)
    l1+frac*(l2-l1)
end
function previous_center(t,nx,lon,lat,mean_single)
    p=t.pixels[end-1]; _,y=xy(p,nx)
    # MATLAB's single input returns a single-precision mean for y.
    my0=t.single[end-1] ? mean(Float32.(y)) : mean(y)
    my=mean_single ? Float32(my0) : Float64(my0)
    mx=mean_lon_index_periodic(p,nx)
    latval=if mean_single
        Float32(lat[floor(Int,my)])+(my-floor(my))*Float32(lat[ceil(Int,my)]-lat[floor(Int,my)])
    else
        lat[floor(Int,my)]+(my-floor(my))*(lat[ceil(Int,my)]-lat[floor(Int,my)])
    end
    latrad=mean_single ? latval*Float32(pi/180) : latval*pi/180
    (interp_lon_periodic(lon,mx)*pi/180,Float64(latrad))
end
function spherical_distance(p,nx,lon,lat,center)
    x=mod1(p,nx); y=fld(p-1,nx)+1
    l1=interp_lon_periodic(lon,x)*pi/180; a1=Float32(lat[y])*Float32(pi/180)
    l2,a2=center; dl=mod(l1-l2+pi,2pi)-pi; da=a1-a2
    # Daily patch coordinates are single in MATLAB. Mixed scalar arithmetic
    # therefore rounds each latitude/alpha operation to Float32, even though
    # longitudes returned by interp_lon_periodic are double.
    da=Float32(a1-Float32(a2))
    a=sin(da/Float32(2))^2+(cos(a1)*Float32(cos(a2)))*Float32(sin(dl/2)^2)
    # abs(sqrt(...)) in MATLAB also handles roundoff outside [0,1].
    angle=Float32(2)*atan(sqrt(abs(a)),sqrt(abs(Float32(1)-a)))
    Float32(6378.137)*angle
end
function patch_components(p,dims)
    mask=falses(dims); mask[p].=true
    bwconncomp_periodic_lon(mask).PixelIdxList
end

function repair_fragments!(search,old,newfusion,dims)
    for (j,id) in enumerate(old)
        p=get(newfusion,j,search[id].pixels[end])
        cc=patch_components(p,dims); length(cc)<=1 && continue
        keep=argmax(length.(cc))
        for k in eachindex(cc)
            k==keep && continue
            fragment=cc[k]
            for (other,target) in enumerate(old)
                other==j && continue
                targetp=search[target].pixels[end]
                joined=vcat(targetp,fragment)
                if length(patch_components(joined,dims))<=length(patch_components(targetp,dims))
                    search[target].pixels[end]=joined; search[target].single[end]=false
                    search[id].pixels[end]=sort!(collect(setdiff(Set(search[id].pixels[end]),Set(fragment))))
                    search[id].single[end]=false
                    break
                end
            end
        end
    end
end

function merge_patches!(search,patches,count,day,dims,lon,lat)
    nx=dims[1]; merged=findall(>(1),count)
    isempty(merged) && return
    oldrows=Vector{Vector{Int}}()
    for patchid in merged
        p=patches[patchid]; xs=mod1.(p,nx); sp=Set(p); ids=Int[]
        for (n,t) in enumerate(search)
            q=t.pixels[end]
            if day==t.day[end] && ((length(p)==length(q) && xs==mod1.(q,nx)) ||
                    sum(intersect(sp,Set(q)))==sum(p))
                push!(ids,n)
            end
        end
        push!(oldrows,ids)
    end
    for (row,patchid) in enumerate(merged)
        old=oldrows[row]; p=patches[patchid]
        isempty(old) && continue
        owner=zeros(Int,length(p)); centers=Tuple{Float64,Float64}[]
        mean_single=search[old[1]].single[end-1]
        for (j,id) in enumerate(old)
            tr=search[id]; overlap=intersect(Set(tr.pixels[end-1]),Set(tr.pixels[end]))
            for k in eachindex(p)
                owner[k]==0 && p[k] in overlap && (owner[k]=j)
            end
            push!(centers,previous_center(tr,nx,lon,lat,mean_single))
        end
        for k in eachindex(p)
            owner[k]!=0 && continue
            owner[k]=argmin([spherical_distance(p[k],nx,lon,lat,c) for c in centers])
        end
        newfusion=Dict{Int,Vector{Int}}()
        for (j,id) in enumerate(old)
            q=search[id].pixels[end]
            if length(q)==length(p)
                # Source uses positional mask on q, even if x-only comparison matched.
                search[id].pixels[end]=q[owner.==j]
            else
                notp=sort!(collect(setdiff(Set(q),Set(p))))
                shared=sort!(collect(intersect(Set(q),Set(p))))
                length(shared)==length(owner) || throw(ArgumentError("MATLAB merge positional mask mismatch"))
                assigned=shared[owner.==j]
                search[id].pixels[end]=vcat(notp,assigned); newfusion[j]=assigned
            end
            length(q)!=length(p) && (search[id].single[end]=false)
        end
        repair_fragments!(search,old,newfusion,dims)
    end
end

"""
    track_mhw_3d(mask, lon, lat, nums=1; min_duration=5, overlap_ratio=0.5)

Faithful daily-patch overlap tracker from the local mhwtrack.m. Each snapshot
uses 8 connectivity and periodic longitude. Splits stay in the parent event;
merges partition the new patch among existing events using prior intersections,
centroid distance, and fragment reassignment. No temporal gaps are bridged.

`nums` is the daily minimum pixel count. Defaults 0.5 and 5 reproduce the fixed
MATLAB constants. Nonzero values including NaN count as active; prepare binary
input explicitly when NaN means no event. Returns Vector{Track3D} in MATLAB
closure order. Empty input returns an empty vector (the source can error).
"""
function track_mhw_3d(mask::AbstractArray{<:Number,3},lon,lat,nums::Real=1;
        min_duration::Real=5,overlap_ratio::Real=0.5)
    nx,ny,nt=size(mask)
    (length(lon),length(lat))==(nx,ny) || throw(DimensionMismatch("coordinate lengths must match mask"))
    search=OpenTrack[]; closed=OpenTrack[]
    for day in 1:nt
        cc=bwconncomp_periodic_lon(@view mask[:,:,day])
        patches=filter(p->length(p)>=nums,cc.PixelIdxList)
        count=zeros(Int,length(patches))
        if day==1
            append!(search,[OpenTrack(day,p,i) for (i,p) in enumerate(patches)])
            continue
        end
        for tr in search
            day-tr.day[end]==1 || continue
            old=Set(tr.pixels[end]); hits=Int[]
            for (i,p) in enumerate(patches)
                denominator=min(length(old),length(p))
                ratio=denominator==0 ? NaN : count_intersection(old,p)/denominator
                ratio>=overlap_ratio && push!(hits,i)
            end
            isempty(hits) && continue
            push!(tr.day,day); push!(tr.pixels,reduce(vcat,patches[hits])); push!(tr.single,true)
            if length(hits)>1
                push!(tr.split_num,length(hits)); push!(tr.split_day,day)
            end
            count[hits].+=1
        end
        merge_patches!(search,patches,count,day,(nx,ny),lon,lat)
        for i in findall(iszero,count)
            # Preserve source ori_order: last patch of identical length.
            order=findlast(p->length(p)==length(patches[i]),patches)
            push!(search,OpenTrack(day,patches[i],order))
        end
        moved=findall(t->t.day[end]<=day-1,search)
        append!(closed,search[moved]); deleteat!(search,moved)
    end
    append!(closed,search)
    [output_track(t,nx) for t in closed if length(t.day)>=min_duration]
end
count_intersection(a,b)=sum(p->p in a,b;init=0)
