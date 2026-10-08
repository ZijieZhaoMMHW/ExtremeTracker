@testset "MATLAB regular 4-D exact reference" begin
    f=matread(joinpath(@__DIR__,"fixtures","lonlat4d.mat"))
    for i in eachindex(f["masks"])
        actual,info=track_mhw_4d_lonlat(f["masks"][i];min_duration=1)
        expected=indexed_from_mat(f["tracks"][i],4)
        test_indexed_exact(actual,expected)
        @test partition_equal(actual,expected)
        for key in (:n_components_before_periodic_union,:n_raw_components,:n_lon_seam_edges,
            :n_active_lon_boundary_cells,:n_lon_seam_component_pairs)
            @test getproperty(info,key)==f["infos"][i][string(key)]
        end
    end
    actual,info=track_mhw_4d_lonlat(f["masks"][1];min_duration=1,time=vec(f["time"]),
        mask_size=Tuple(Int.(vec(f["global_size"]))),add_idx5_alias=true)
    test_indexed_exact(actual,indexed_from_mat(f["mapped"],4))
    @test all(t->t.idx5===t.idx4,actual)
end
@testset "MATLAB LLC 4-D exact reference" begin
    f=matread(joinpath(@__DIR__,"fixtures","llc4d.mat"))
    topology=LLCTopology(f["seamEdges"],(3,2,2))
    for i in eachindex(f["masks"])
        actual,info=track_mhw_4d_llc(f["masks"][i],topology;min_duration=1)
        expected=indexed_from_mat(f["tracks"][i],5)
        test_indexed_exact(actual,expected)
        @test partition_equal(actual,expected)
        for key in (:n_tile_components,:n_raw_components,:n_active_seam_cells,:n_seam_component_pairs)
            @test getproperty(info,key)==f["infos"][i][string(key)]
        end
    end
end
@testset "4-D edge and filtering rules" begin
    a=falses(5,4,3,5); a[1,1,1,1]=a[5,2,2,2]=true
    @test length(track_mhw_4d_lonlat(a;min_duration=1)[1])==1
    a[5,2,2,2]=false; a[5,2,3,2]=true
    @test length(track_mhw_4d_lonlat(a;min_duration=1)[1])==2
    @test isempty(track_mhw_4d_lonlat(falses(2,2,1,5))[1])
    @test length(track_mhw_4d_lonlat(trues(2,2,2,5))[1])==1
    a=ones(2,2,1,5)
    @test isempty(track_mhw_4d_lonlat(a;time=fill(1,5))[1])
    @test_throws ArgumentError track_mhw_4d_lonlat(a;time=[0,1,2,3,4])
    @test_throws ArgumentError track_mhw_4d_lonlat(ones(1,2,2,5))
    @test_throws DimensionMismatch track_mhw_4d_lonlat(a;mask_size=(3,2,1,5))
    # Tile numbers have no implicit array-axis adjacency.
    a=trues(2,2,2,1,2)
    @test length(track_mhw_4d_llc(a,zeros(Int,0,2);min_duration=1)[1])==2
    @test_throws ArgumentError LLCTopology([1 1],(2,2,2))
    @test_throws ArgumentError LLCTopology([1 5;5 1],(2,2,2))
    @test_throws ArgumentError LLCTopology([1 2],(2,2,2))
    @test_throws ArgumentError LLCTopology([0 5],(2,2,2))
    # Nonfinite numeric input follows ~= 0, no invented land rule.
    a=zeros(2,2,1,1); a[1]=NaN
    @test length(track_mhw_4d_lonlat(a;min_duration=1)[1])==1
end
