"""
    load_test_data(directory; binary_rule=:raw)

Read the inspected mhw_ts.mat and lonlatoi.mat without permuting dimensions.
`:raw` preserves all values including NaNs; `:positive` explicitly maps finite
positive values to true; `:finite` maps every finite value to true.
"""
function load_test_data(directory;binary_rule=:raw)
    a=matread(joinpath(directory,"mhw_ts.mat"))["mhw_ts"]
    coords=matread(joinpath(directory,"lonlatoi.mat"))
    ndims(a)==3 || throw(DimensionMismatch("mhw_ts must have lon×lat×time dimensions"))
    lon=vec(coords["lon_used"]); lat=vec(coords["lat_used"])
    size(a)[1:2]==(length(lon),length(lat)) || throw(DimensionMismatch("coordinates do not match data"))
    mask=if binary_rule==:raw; a
    elseif binary_rule==:positive; isfinite.(a).&(a.>0)
    elseif binary_rule==:finite; isfinite.(a)
    else; throw(ArgumentError("binary_rule must be :raw, :positive, or :finite"))
    end
    (mask=mask,lon=lon,lat=lat,binary_rule=binary_rule)
end
"""
    load_matlab_tracks(path; variable="tracks")

Read a MATLAB scalar struct or struct array into `Vector{Track3D}`. Cell arrays
become one coordinate vector per day. Preserve all seven track fields, one-based
indices, column-major order, and absolute day labels; no transpose or reindexing.
"""
function load_matlab_tracks(path; variable="tracks")
    raw=matread(path)[variable]
    records=raw isa AbstractDict ? [raw] : raw isa MAT.MatlabStructArray ?
        vec(Array{Dict{String,Any}}(raw)) : vec(raw)
    asvec(x::Number)=[x]
    asvec(x)=vec(x)
    result=Track3D[]
    for s in records
        days=Int.(asvec(s["day"]))
        xs=[Float64.(asvec(p)) for p in asvec(s["xloc"])]
        ys=[Float64.(asvec(p)) for p in asvec(s["yloc"])]
        length(days)==length(xs)==length(ys) || throw(DimensionMismatch("track day/cell counts differ"))
        all(length(x)==length(y) for (x,y) in zip(xs,ys)) || throw(DimensionMismatch("daily xloc/yloc lengths differ"))
        push!(result,Track3D(days,xs,ys,Int(s["ori_day"]),Int(s["ori_order"]),
            Int.(asvec(s["split_num"])),Int.(asvec(s["split_day"]))))
    end
    result
end
function matlab_struct(t::Track3D)
    Dict("day"=>reshape(t.day,:,1),"xloc"=>reshape(Any[t.xloc...],:,1),
         "yloc"=>reshape(Any[t.yloc...],:,1),"ori_day"=>t.ori_day,"ori_order"=>t.ori_order,
         "split_num"=>reshape(t.split_num,:,1),"split_day"=>reshape(t.split_day,:,1))
end
function matlab_struct(t::IndexedTrack{N}) where N
    d=Dict{String,Any}("day"=>reshape(t.day,:,1),"idx$N"=>reshape(Any[t.indices...],:,1),"nvoxels"=>reshape(t.nvoxels,:,1))
    t.has_idx5_alias && (d["idx5"]=d["idx$N"])
    d
end
function save_tracks(path,tracks,info)
    payload=if isempty(tracks)
        T=eltype(tracks)
        names=if T==Track3D
            collect(string.(fieldnames(Track3D)))
        else
            N=T.parameters[1]
            fields=["day","idx$N","nvoxels"]
            hasproperty(info,:has_idx5_alias) && info.has_idx5_alias && push!(fields,"idx5")
            fields
        end
        MAT.MatlabStructArray(names,(0,1))
    else
        reshape([matlab_struct(t) for t in tracks],:,1)
    end
    matwrite(path,Dict("tracks"=>payload,
        "info"=>Dict(string(k)=>(v isa Tuple ? collect(v) : v) for (k,v) in pairs(info))))
end
