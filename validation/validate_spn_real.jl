using MHWTracking, MAT, Test
include(joinpath(@__DIR__,"..","test","reference_helpers.jl"))
@testset "Actual field spn MATLAB equivalence" begin
    for flag in 0:1
        f=matread(joinpath(@__DIR__,"..","test","fixtures","spn_real_$flag.mat"))
        r,s,n,x,y=spn(overlap_from_mat(f["tr"]),vec(f["lon"]),vec(f["lat"]),f["anomaly"],flag)
        for i in eachindex(s)
            a=s[i]; b=f["snapshots"][i]; finite=isfinite.(a).&isfinite.(b)
            println("flag=",flag," event=",i," NaN mismatches=",count(isnan.(a).!=isnan.(b)),
                " zero mismatches=",count((a.==0).!=(b.==0))," max finite error=",
                any(finite) ? maximum(abs.(a[finite].-b[finite])) : "no common finite values")
            mismatch=findall(isnan.(a).!=isnan.(b))
            for p in mismatch[1:min(end,45)]
                println(p," Julia=",a[p]," MATLAB=",b[p])
            end
            @test float_match(a,b)
        end
        @test float_match(n,reshape(f["normalized"],size(n)))
    end
end
