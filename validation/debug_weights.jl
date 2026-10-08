using MHWTracking, MAT, Test
import DelaunayTriangulation as DT
include(joinpath(@__DIR__,"..","test","reference_helpers.jl"))
f=matread(joinpath(@__DIR__,"..","test","fixtures","spn_1.mat"))
tr=overlap_from_mat(f["tr"])[1]; lon=vec(f["lon"]); lat=vec(f["lat"])
xs=Int.(tr.xloc[1]); ys=Int.(tr.yloc[1]); cx=sum(lon[xs])/length(xs);cy=sum(lat[ys])/length(ys)
r=maximum(MHWTracking.geodist(lat[ys[k]],cy,lon[xs[k]],cx)[1] for k in eachindex(xs))
# The event-wide maximum can be on another day; use the MATLAB radius.
r=f["radius"]; x,y=MHWTracking.polar_grid(r)
allpoints=[MHWTracking.geodist(cy,yy,cx,xx)[2:3] for yy in lat for xx in lon]
circle=[(x[100,j],y[100,j]) for j in 1:100]
ids=findall(q->MHWTracking.inpolygon(q,circle),allpoints); ps=allpoints[ids]
vals=vec(f["anomaly"][:,:,1])[ids]; tri=DT.triangulate(ps;randomise=false)
println("center ",(cx,cy))
for row in 1:6,col in (1,100)
 q=(x[row,col],y[row,col]); a,b,c=DT.find_triangle(tri,q;k=1)
 den=MHWTracking.cross2(ps[a],ps[b],ps[c])
 w=(MHWTracking.cross2(q,ps[b],ps[c])/den,MHWTracking.cross2(ps[a],q,ps[c])/den,
    MHWTracking.cross2(ps[a],ps[b],q)/den)
 println((row,col)," q=",q," triangle=",(ps[a],ps[b],ps[c])," values=",(vals[a],vals[b],vals[c])," w=",w," matlab=",f["snapshots"][1][row,col,1])
end
