@testset "Neighborhood and smoothing MATLAB reference" begin
    f=matread(joinpath(@__DIR__,"fixtures","auxiliary.mat"))
    near=grid_nearby_index(vec(f["lon_aux"]),vec(f["lat_aux"]),f["res"])
    @test all(near[k]==vecmat(f["near"][k]) for k in eachindex(near))
    actual=convzz(f["data"],near,f["ocean"])
    @test float_match(actual,f["smooth"])
    @test isnan(actual[3,2])
    @test all(isnan,convzz(fill(NaN,6,3),near,trues(6,3)))
    @test_throws ArgumentError grid_nearby_index([0,180],[-1,1],400)
end
@testset "spn real MATLAB numeric fixtures" begin
    for flag in 0:1
        f=matread(joinpath(@__DIR__,"fixtures","spn_$flag.mat"))
        tracks=overlap_from_mat(f["tr"])
        r,s,n,x,y=spn(tracks,vec(f["lon"]),vec(f["lat"]),f["anomaly"],flag)
        @test isapprox(r,vecmat(f["radius"]);atol=1e-10,rtol=1e-12)
        @test float_match(x,f["xunit"];atol=2e-15)
        @test float_match(y,f["yunit"];atol=2e-15)
        @test float_match(s[1],f["snapshots"][1])
        # MATLAB drops the singleton fifth dimension in v5 MAT storage.
        @test float_match(n,reshape(f["normalized"],size(n)))
    end
end
@testset "spn actual data coordinates, tracks and field" begin
    for flag in 0:1
        f=matread(joinpath(@__DIR__,"fixtures","spn_real_$flag.mat"))
        tracks=overlap_from_mat(f["tr"])
        r,s,n,x,y=spn(tracks,vec(f["lon"]),vec(f["lat"]),f["anomaly"],flag)
        @test isapprox(r,vecmat(f["radius"]);atol=1e-10,rtol=1e-12)
        for i in eachindex(s)
            @test float_match(s[i],f["snapshots"][i])
        end
        @test float_match(n,reshape(f["normalized"],size(n)))
    end
end
@testset "Indexed polygon membership preserves the original predicate" begin
    rng=MersenneTwister(6150)
    polygons=[Tuple{Float64,Float64}[],[(0.,0.)],[(0.,0.),(2.,0.)],
        [(0.,0.),(3.,0.),(3.,1.),(1.,1.),(1.,3.),(0.,3.)],
        [(0.,0.),(1.,0.),(1.,1.),(0.,1.),(0.,0.),(-1.,0.),(-1.,-1.),(0.,-1.)],
        [(3433sin(t),3433cos(t)) for t in range(0,2pi;length=100)]]
    for polygon in polygons
        lookup=MHWTracking.PolygonLookup(polygon)
        bound=maximum((max(abs(p[1]),abs(p[2])) for p in polygon);init=1.)
        queries=[(4bound*rand(rng)-2bound,4bound*rand(rng)-2bound) for _ in 1:1000]
        append!(queries,polygon)
        for i in eachindex(polygon)
            a=polygon[i]; b=polygon[mod1(i+1,length(polygon))]
            midpoint=((a[1]+b[1])/2,(a[2]+b[2])/2)
            push!(queries,midpoint,(midpoint[1]+1e-10,midpoint[2]-1e-10))
        end
        @test all(MHWTracking.inpolygon(q,lookup)==MHWTracking.inpolygon(q,polygon) for q in queries)
    end
end
@testset "Tasman Sea holes touching the exterior at a vertex" begin
    f=matread(joinpath(@__DIR__,"fixtures","boundary_tasman.mat"))
    for geometry in vec(f["geometries"])
        points=[(p[1],p[2]) for p in eachrow(geometry["points"])]
        boundary=MHWTracking.boundary_compact(points)
        lookup=MHWTracking.PolygonLookup(boundary)
        actual=map((x,y)->MHWTracking.inpolygon((x,y),lookup),f["x"],f["y"])
        @test actual==geometry["inside"]
        @test all(p->MHWTracking.inpolygon(p,lookup),points)
    end
end
