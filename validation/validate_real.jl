using MHWTracking, MAT, Test, Serialization
include(joinpath(@__DIR__,"..","test","reference_helpers.jl"))
directory=length(ARGS)>0 ? ARGS[1] : get(ENV,"MHW_DATA_DIR",normpath(joinpath(@__DIR__,"..")))
reference=length(ARGS)>1 ? ARGS[2] : joinpath(@__DIR__,"results","matlab_real_positive.mat")
data=load_test_data(directory;binary_rule=:positive)
println("Loaded positive binary mask: ",size(data.mask),"; active=",count(data.mask)); flush(stdout)
cache=joinpath(@__DIR__,"results","julia_real_cache.bin")
if "--cached" in ARGS
    saved=deserialize(cache)
    actual=saved.tracks; seconds=saved.seconds; allocated=saved.allocated_bytes
else
    # Compile on an independent small input before measuring the full workload.
    track_mhw_3d(data.mask[1:10,1:10,1:6],data.lon[1:10],data.lat[1:10],4)
    GC.gc()
    measurement=@timed track_mhw_3d(data.mask,data.lon,data.lat,4)
    actual=measurement.value; seconds=measurement.time; allocated=measurement.bytes
    serialize(cache,(tracks=actual,seconds=seconds,allocated_bytes=allocated))
end
println("Julia: seconds=",seconds," allocated_bytes=",allocated," tracks=",length(actual)); flush(stdout)
if !isfile(reference)
    println("MATLAB reference still pending; Julia result cached for comparison.")
    exit()
end
f=matread(reference); expected=overlap_from_mat(f["tracks"])
exact=length(actual)==length(expected) && all(all(getfield(a,k)==getfield(b,k) for k in fieldnames(Track3D)) for (a,b) in zip(actual,expected))
partition=partition_equal(actual,expected,size(data.mask))
result=Dict("shape"=>collect(size(data.mask)),"active"=>count(data.mask),"nums"=>4,
    "tracks_julia"=>length(actual),"tracks_matlab"=>length(expected),"strict_fields_equal"=>exact,
    "event_partition_equal"=>partition,"seconds"=>seconds,"allocated_bytes"=>allocated,
    "matlab_seconds"=>f["s"]["elapsed_seconds"],"julia_version"=>string(VERSION))
matwrite(joinpath(@__DIR__,"results","julia_real_metrics.mat"),result)
open(joinpath(@__DIR__,"results","real_metrics.txt"),"w") do io
    for k in sort!(collect(keys(result))); println(io,k," = ",result[k]); end
end
println(result)
@test exact
@test partition
