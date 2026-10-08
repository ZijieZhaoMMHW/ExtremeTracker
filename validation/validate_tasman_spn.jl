using MHWTracking, MAT, Test
include(joinpath(@__DIR__,"..","test","reference_helpers.jl"))
include(joinpath(@__DIR__,"..","examples","tasman_data.jl"))
root=normpath(joinpath(@__DIR__,".."))
directory=get(ENV,"MHW_TASMAN_DATA_DIR",joinpath(root,"data","tasman_sea"))
data=load_tasman_data(directory)
report=String[]
@testset "Full 230-day Tasman Sea MATLAB equivalence" begin
    @test length(data.track.day)==230
    @test data.track.day==collect(12329:12558)
    @test length(data.track.split_day)==38
    for flag in 0:1
        measured=run_or_load_tasman_spn(data,flag;output_dir=joinpath(@__DIR__,"results"),recompute="--recompute" in ARGS)
        actual=measured.output
        f=matread(joinpath(@__DIR__,"results","matlab_tasman_spn_$flag.mat"))
        expected=(vecmat(f["radius"]),[only(vec(f["snapshots"]))],reshape(f["normalized"],100,100,5,1,1),f["xunit"],f["yunit"])
        @test isapprox(actual[1],expected[1];atol=1e-10,rtol=1e-12)
        max_error=0.0; nan_mismatch=0; zero_mismatch=0
        for (a,b) in ((only(actual[2]),only(expected[2])),(actual[3],expected[3]),(actual[4],expected[4]),(actual[5],expected[5]))
            @test size(a)==size(b)
            nan_mismatch+=count(isnan.(a).!=isnan.(b))
            zero_mismatch+=count((a.==0).!=(b.==0))
            both=isfinite.(a).&isfinite.(b)
            max_error=max(max_error,any(both) ? maximum(abs.(a[both].-b[both])) : 0.0)
            @test float_match(a,b)
        end
        @test nan_mismatch==0
        @test zero_mismatch==0
        line="flag=$flag radius_km=$(only(actual[1])) max_abs_error=$max_error nan_mismatch=$nan_mismatch zero_mismatch=$zero_mismatch julia_seconds=$(measured.seconds) matlab_seconds=$(f["elapsed_seconds"])"
        println(line); push!(report,line)
    end
    f=matread(joinpath(@__DIR__,"results","matlab_tasman_smoothing.mat"))
    near=grid_nearby_index(vec(f["lon"]),vec(f["lat"]),2.5)
    ocean=matread(joinpath(directory,"ocean_idx.mat"))["ocean_idx"]
    binary=Float64.(.!isnan.(matread(joinpath(directory,"snaps.mat"))["snaps"]))
    binary[.!ocean].=NaN
    fraction=convzz(binary,near,ocean)
    @test all(near[k]==vecmat(f["near"][k]) for k in eachindex(near))
    @test float_match(binary,f["binary"])
    @test float_match(fraction,f["fraction"])
    @test Float64.(fraction.>=0.5)==f["smoothed"]
end
write(joinpath(@__DIR__,"results","tasman_metrics.txt"),join(report,"\n")*"\n")
