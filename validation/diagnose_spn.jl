using MAT, MHWTracking, Test
include(joinpath(@__DIR__,"..","test","reference_helpers.jl"))
for flag in 0:1
    f=matread(joinpath(@__DIR__,"..","test","fixtures","spn_$flag.mat"))
    r,s,n,x,y=spn(overlap_from_mat(f["tr"]),vec(f["lon"]),vec(f["lat"]),f["anomaly"],flag)
    b=f["snapshots"][1]; a=s[1]; finite=isfinite.(a).&isfinite.(b)
    println("flag=",flag," radius=",r," MATLAB=",f["radius"])
    println("NaN mismatch=",count(isnan.(a).!=isnan.(b))," finite maxerr=",maximum(abs.(a[finite].-b[finite])))
    println("Zero-mask mismatch=",count((a.==0).!=(b.==0)))
    for p in findall(isnan.(a).!=isnan.(b))[1:min(end,8)]
        println(p," julia=",a[p]," matlab=",b[p])
    end
    matwrite(joinpath(@__DIR__,"results","spn_julia_$flag.mat"),Dict("snapshots"=>a,"normalized"=>n))
end
