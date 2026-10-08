nanmean(v) = begin
    s=0.0; n=0
    for x in v
        if !isnan(x); s+=x; n+=1; end
    end
    n==0 ? NaN : s/n
end

"""
    grid_nearby_index(lon, lat, res)

Coordinate window ±res, with the MATLAB 0..360 longitude boundary logic and
bounded latitude. Cell contents retain meshgrid(:) ordering. A window crossing
both ends is undefined in MATLAB and is rejected here.
"""
function grid_nearby_index(lon,lat,res::Real)
    nx=length(lon); ny=length(lat); out=Matrix{Vector{Int}}(undef,nx,ny)
    for x in 1:nx
        lo=lon[x]-res; hi=lon[x]+res
        pred = if lo>=0 && hi<360
            v->lo<=v<=hi
        elseif lo<0 && hi<360
            v->(0<=v<=hi || 360+lo<=v<360)
        elseif lo>=0 && hi>=360
            v->(lo<=v<360 || 0<=v<=hi-360)
        else
            throw(ArgumentError("window crosses both longitude bounds; undefined in MATLAB"))
        end
        xs=findall(pred,vec(lon))
        for y in 1:ny
            ys=findall(v->lat[y]-res<=v<=lat[y]+res,vec(lat))
            out[x,y]=[xx+(yy-1)*nx for xx in xs for yy in ys]
        end
    end
    out
end

"""
    convzz(data, indices, ocean_mask)

NaN-ignoring unweighted neighborhood mean at nonzero ocean centers; other
centers remain NaN. Neighbor values are not additionally masked by ocean_mask.
"""
function convzz(data::AbstractMatrix,indices,ocean_mask)
    size(data)==size(indices)==size(ocean_mask) || throw(DimensionMismatch("grid shapes must match"))
    out=fill(NaN,size(data))
    for p in eachindex(data)
        active(ocean_mask[p]) || continue
        out[p]=nanmean(data[q] for q in indices[p])
    end
    out
end
