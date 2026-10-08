using MAT,MHWTracking,Test
import DelaunayTriangulation as DT
f=matread(joinpath(@__DIR__,"..","test","fixtures","boundary_audit.mat"))
for i in axes(f["geometries"],1),j in axes(f["geometries"],2)
    g=f["geometries"][i,j]; g isa Dict || continue
    points=[Tuple(row) for row in eachrow(g["points"])]
    p=MHWTracking.boundary_compact(points)
    expected=points[Int.(vec(g["indices"]))]
    println("event=",i," day=",g["day"]," points=",length(points)," MATLAB alpha=",f["critical"][i,j],
        " max=",f["maxalpha"][i,j]," Julia border=",length(p)," MATLAB border=",length(expected),
        " vertex set difference=",length(symdiff(Set(p),Set(expected))))
    if Set(p)!=Set(expected)
        println("Julia extra indices=",findall(q->q in setdiff(Set(p),Set(expected)),points))
        println("MATLAB extra indices=",findall(q->q in setdiff(Set(expected),Set(p)),points))
    end
    function edge_set(poly)
        Set(minmax(findfirst(==(poly[k]),points),findfirst(==(poly[mod1(k+1,length(poly))]),points))
            for k in eachindex(poly) if poly[k]!=poly[mod1(k+1,length(poly))])
    end
    if edge_set(p)!=edge_set(expected)
        println("Julia extra edges=",setdiff(edge_set(p),edge_set(expected)))
        println("MATLAB extra edges=",setdiff(edge_set(expected),edge_set(p)))
        ts=sort!(collect(DT.each_solid_triangle(DT.triangulate(points;randomise=false))))
        rs=[MHWTracking.circumradius(points[a],points[b],points[c]) for (a,b,c) in ts]
        ds=MHWTracking.DSU(length(ts)); used=falses(length(points)); accepted=falses(length(ts)); es=Dict{Tuple{Int,Int},Int}()
        for k in sortperm(rs)
            accepted[k]=true; a,b,c=ts[k]; used[[a,b,c]].=true
            for e in (minmax(a,b),minmax(b,c),minmax(c,a))
                if haskey(es,e); MHWTracking.union!(ds,k,es[e]); else; es[e]=k; end
            end
            roots=Set(MHWTracking.root!(ds,l) for l in eachindex(accepted) if accepted[l])
            if all(used) && length(roots)==1
                println("Julia computed critical=",rs[k]," triangle=",ts[k]); break
            end
        end
        println("Near-critical radii=",sort!([(rs[k],ts[k]) for k in eachindex(rs) if abs(rs[k]-f["critical"][i,j])<1e-6]))
    end
end
