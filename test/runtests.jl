using Test, Random, MAT, MHWTracking
include("reference_helpers.jl")
@testset "MHWTracking" begin
    include("test_tracking_3d.jl")
    include("test_tracking_4d.jl")
    include("test_auxiliary.jl")
    include("test_seam_edges.jl")
    include("test_io.jl")
end
