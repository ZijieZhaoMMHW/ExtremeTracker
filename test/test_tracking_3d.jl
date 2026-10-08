@testset "MATLAB overlap reference" begin
    f=matread(joinpath(@__DIR__,"fixtures","overlap.mat"))
    for i in eachindex(f["inputs"])
        isempty(f["errors"][i]) || continue # recorded source failures are not success claims
        @testset "$(f["names"][i])" begin
            mask=f["inputs"][i]
            actual=track_mhw_3d(mask,vec(f["lon"]),vec(f["lat"]))
            expected=overlap_from_mat(f["outputs"][i])
            test_overlap_exact(actual,expected)
            @test partition_equal(actual,expected,size(mask))
        end
    end
end
@testset "3-D boundaries and method distinction" begin
    @test isempty(track_mhw_3d(falses(4,3,6),1:4,1:3))
    @test isempty(track_mhw_3d(falses(4,3,0),1:4,1:3))
    a=ones(4,3,6); t=track_mhw_3d(a,1:4,1:3)
    @test length(t)==1 && length(t[1].day)==6
    @test all(labels(t,size(a)).==1)
    # No overlap, but diagonal space-time connection.
    a=falses(4,3,2); a[2,2,1]=a[3,2,2]=true
    @test length(track_mhw_3d(a,1:4,1:3;min_duration=1))==2
    @test length(track_mhw_3d_connected(a;min_duration=1)[1])==1
    a=falses(5,4,2); a[1,1,1]=a[5,2,2]=true
    @test bwconncomp_periodic_lon(a).NumObjects==1
    @test bwconncomp_periodic_lon(a;connectivity=6).NumObjects==2
    a=zeros(5,4,5); a[2,1,:].=NaN
    @test length(track_mhw_3d(a,1:5,1:4))==1
    @test bwconncomp_periodic_lon(ones(1,1)).NumObjects==1
    # Periodic helper concatenates pre-union components, not globally sorted pixels.
    a=falses(5,3); a[1,2]=a[5,1]=true
    @test bwconncomp_periodic_lon(a).PixelIdxList==[[5,6]]
    @test_throws DimensionMismatch track_mhw_3d(ones(4,3,6),1:3,1:3)
end
