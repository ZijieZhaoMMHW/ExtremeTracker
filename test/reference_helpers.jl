vecmat(x::Number)=[x]
vecmat(x)=vec(x)
function structvec(x)
    x isa Dict && return [x]
    isempty(x) && return Dict{String,Any}[]
    x isa MAT.MatlabStructArray && return vec(Array{Dict{String,Any}}(x))
    vec(x)
end
function overlap_from_mat(x)
    [Track3D(Int.(vecmat(s["day"])),[Float64.(vecmat(p)) for p in vecmat(s["xloc"])],
      [Float64.(vecmat(p)) for p in vecmat(s["yloc"])],Int(get(s,"ori_day",0)),Int(get(s,"ori_order",0)),
      Int.(vecmat(get(s,"split_num",Int[]))),Int.(vecmat(get(s,"split_day",Int[])))) for s in structvec(x)]
end
function indexed_from_mat(x,N)
    [IndexedTrack{N}(Int32.(vecmat(s["day"])),[UInt64.(vecmat(p)) for p in vecmat(s["idx$N"])],
        UInt32.(vecmat(s["nvoxels"])),haskey(s,"idx5") && N!=5) for s in structvec(x)]
end
function test_overlap_exact(a,b)
    @test length(a)==length(b)
    length(a)==length(b) || return
    for (x,y) in zip(a,b), key in fieldnames(Track3D)
        @test getfield(x,key)==getfield(y,key)
    end
end
function test_indexed_exact(a,b)
    @test length(a)==length(b)
    length(a)==length(b) || return
    for (x,y) in zip(a,b)
        @test x.day==y.day
        @test x.indices==y.indices
        @test x.nvoxels==y.nvoxels
    end
end
function float_match(a,b;atol=1e-10,rtol=1e-11)
    size(a)==size(b) && isnan.(a)==isnan.(b) &&
        all(isapprox.(a[.!isnan.(a)],b[.!isnan.(b)];atol,rtol))
end
