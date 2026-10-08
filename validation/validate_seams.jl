using MHWTracking, MAT, Test
path=length(ARGS)>0 ? ARGS[1] : joinpath(@__DIR__,"..","test","fixtures","llc90_topology.mat")
topology=load_seam_edges(path); edges=topology.seam_edges; grid=topology.grid_size
n2=prod(grid); nz=3; batch_size=16
cases=(("same",1,1,0,0,1),("forwardDiagonal",1,2,0,1,1),
    ("reverseDiagonal",2,1,1,0,1),("depthGapTwo",1,3,0,0,2),("timeGapTwo",1,1,0,2,2))
start=time()
@testset "All $(size(edges,1)) actual LLC90 seam edges" begin
    for (name,za,zb,ta,tb,expected) in cases
        @testset "$name" begin
            for first in 1:batch_size:size(edges,1)
                last=min(first+batch_size-1,size(edges,1)); n=last-first+1
                dims=(grid...,nz,4n); ps=Int[]; expected_memberships=Vector{Vector{UInt64}}()
                for (j,row) in enumerate(first:last)
                    t0=1+4(j-1); a,b=edges[row,:]
                    ia=a+(za-1)*n2+(t0+ta-1)*n2*nz
                    ib=b+(zb-1)*n2+(t0+tb-1)*n2*nz
                    append!(ps,[ia,ib])
                    if expected==1; push!(expected_memberships,sort!(UInt64[ia,ib]))
                    else; push!(expected_memberships,[UInt64(ia)]); push!(expected_memberships,[UInt64(ib)])
                    end
                end
                tracks,info=track_mhw_4d_llc(SparseMask(dims,ps),topology;min_duration=1)
                @test info.n_raw_components==expected*n
                @test event_partition(tracks)==sort!(expected_memberships)
            end
        end
        println(name," passed for ",size(edges,1)," edges"); flush(stdout)
    end
end
open(joinpath(@__DIR__,"results","seam_metrics.txt"),"w") do io
    println(io,"edges = ",size(edges,1)); println(io,"cases = ",length(cases))
    println(io,"isolated_pairs = ",size(edges,1)*length(cases)); println(io,"seconds = ",time()-start)
end
