@testset "Portable exact LLC seam cases" begin
    topo=LLCTopology([3 10;6 7],(3,2,2))
    ci=CartesianIndices(topo.grid_size)
    for pair in eachrow(topo.seam_edges), (za,zb,ta,tb,expected) in
        ((1,1,1,1,1),(1,2,1,2,1),(2,1,2,1,1),(1,3,1,1,2),(1,1,1,3,2))
        mask=falses(3,2,2,3,4)
        mask[Tuple(ci[pair[1]])...,za,ta]=true
        mask[Tuple(ci[pair[2]])...,zb,tb]=true
        tracks,info=track_mhw_4d_llc(mask,topo;min_duration=1)
        @test info.n_raw_components==expected
        @test sum(sum(t.nvoxels) for t in tracks)==2
    end
end
@testset "Actual LLC90 topology and grid" begin
    f=matread(joinpath(@__DIR__,"fixtures","llc90_topology.mat"))
    topo=LLCTopology(f["seamEdges"])
    @test size(topo.seam_edges)==(2312,2)
    @test size(f["XC"])==size(f["YC"])==(90,90,13)
    @test length(f["Z"])==50
    @test all(isfinite,f["XC"]) && all(isfinite,f["YC"])
    @test all(topo.seam_edges.==f["seamEdges"])
    # Cover first/last rows and actual corner cells with the public tracker.
    ci=CartesianIndices(topo.grid_size)
    corners=findall(r->any(p->ci[p][1] in (1,90) && ci[p][2] in (1,90),topo.seam_edges[r,:]),1:2312)
    for row in unique(vcat([1,2312],corners))
        a,b=topo.seam_edges[row,:]; n2=prod(topo.grid_size)
        mask=SparseMask((90,90,13,2,2),[a,b+n2+n2*2])
        tracks,_=track_mhw_4d_llc(mask,topo;min_duration=1)
        @test event_partition(tracks)==[sort!(UInt64[a,b+n2+n2*2])]
    end
end
