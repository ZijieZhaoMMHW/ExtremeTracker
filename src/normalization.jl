# Equirectangular distance, exactly the local geodist formula; no lon wrap.
function geodist(lat1,lat2,lon1,lon2)
    a=lat1*pi/180; b=lat2*pi/180
    x=6371*(lon2*pi/180-lon1*pi/180)*cos((a+b)/2)
    y=6371*(b-a)
    hypot(x,y),x,y
end
function polar_grid(radius)
    r=collect(range(0.0,radius;length=100)); theta=collect(range(0.0,2pi;length=100))
    (r.*transpose(sin.(theta)),r.*transpose(cos.(theta)))
end
cross2(a,b,c)=(b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1])
function inpolygon(q,polygon)
    length(polygon)<2 && return false
    inside=false; x,y=q
    for i in eachindex(polygon)
        a=polygon[i]; b=polygon[mod1(i+1,length(polygon))]
        tol=1e-11*max(1.0,abs(x),abs(y),abs(a[1]),abs(a[2]),abs(b[1]),abs(b[2]))
        if abs(cross2(a,b,q))<=tol*max(1.0,hypot(b[1]-a[1],b[2]-a[2])) &&
            min(a[1],b[1])-tol<=x<=max(a[1],b[1])+tol && min(a[2],b[2])-tol<=y<=max(a[2],b[2])+tol
            return true
        end
        if (a[2]>y)!=(b[2]>y)
            x<(b[1]-a[1])*(y-a[2])/(b[2]-a[2])+a[1] && (inside=!inside)
        end
    end
    inside
end
function circumradius(a,b,c)
    area2=abs(cross2(a,b,c)); area2==0 && return Inf
    hypot(a[1]-b[1],a[2]-b[2])*hypot(a[1]-c[1],a[2]-c[2])*hypot(b[1]-c[1],b[2]-c[2])/(2area2)
end

# Index the original polygon edges by their y intervals. The membership formula
# and point-specific tolerance below are identical to inpolygon; the index only
# excludes edges that cannot cross the ray or contain the query point.
struct PolygonLookup
    points::Vector{Tuple{Float64,Float64}}
    bins::Vector{Vector{Int}}
    bounds::NTuple{4,Float64}
    margin::Float64
    scale::Float64
end
function PolygonLookup(polygon)
    points=Tuple{Float64,Float64}.(polygon)
    isempty(points) && return PolygonLookup(points,[Int[]],(0.,0.,0.,0.),0.,1.)
    xmin,xmax=extrema(first.(points)); ymin,ymax=extrema(last.(points))
    margin=4e-11*max(1.,abs(xmin),abs(xmax),abs(ymin),abs(ymax))
    n=clamp(length(points),16,512)
    scale=n/(ymax-ymin+2margin)
    bins=[Int[] for _ in 1:n]
    for i in eachindex(points)
        ya=points[i][2]; yb=points[mod1(i+1,length(points))][2]
        lo=clamp(floor(Int,(min(ya,yb)-ymin)*scale)+1,1,n)
        hi=clamp(floor(Int,(max(ya,yb)-ymin+2margin)*scale)+1,1,n)
        for bin in lo:hi; push!(bins[bin],i); end
    end
    PolygonLookup(points,bins,(xmin,xmax,ymin,ymax),margin,scale)
end
function inpolygon(q,lookup::PolygonLookup)
    points=lookup.points; length(points)<2 && return false
    x,y=q; xmin,xmax,ymin,ymax=lookup.bounds; margin=lookup.margin
    xmin-margin<=x<=xmax+margin && ymin-margin<=y<=ymax+margin || return false
    bin=clamp(floor(Int,(y-ymin+margin)*lookup.scale)+1,1,length(lookup.bins))
    inside=false
    for i in lookup.bins[bin]
        a=points[i]; b=points[mod1(i+1,length(points))]
        tol=1e-11*max(1.0,abs(x),abs(y),abs(a[1]),abs(a[2]),abs(b[1]),abs(b[2]))
        if abs(cross2(a,b,q))<=tol*max(1.0,hypot(b[1]-a[1],b[2]-a[2])) &&
            min(a[1],b[1])-tol<=x<=max(a[1],b[1])+tol && min(a[2],b[2])-tol<=y<=max(a[2],b[2])+tol
            return true
        end
        if (a[2]>y)!=(b[2]>y)
            x<(b[1]-a[1])*(y-a[2])/(b[2]-a[2])+a[1] && (inside=!inside)
        end
    end
    inside
end
# MATLAB HoleThreshold fills every bounded complement component. Sharing
# only a vertex with the exterior does not open a hole: face connectivity
# across shared edges distinguishes those contacts without numeric tolerances.
function alpha_exterior(triangles,radii,alpha)
    edgefaces=Dict{Tuple{Int,Int},Vector{Int}}()
    for (k,(a,b,c)) in enumerate(triangles), e in ((a,b),(b,c),(c,a))
        push!(get!(edgefaces,minmax(e...),Int[]),k)
    end
    selected=radii.<=alpha
    adjacent=[Int[] for _ in triangles]; outside=falses(length(triangles)); queue=Int[]
    for faces in values(edgefaces)
        if length(faces)==1
            k=only(faces)
            if !selected[k] && !outside[k]; outside[k]=true;push!(queue,k);end
        else
            a,b=faces
            push!(adjacent[a],b);push!(adjacent[b],a)
        end
    end
    head=1
    while head<=length(queue)
        k=queue[head];head+=1
        for n in adjacent[k]
            if !selected[n] && !outside[n]
                outside[n]=true;push!(queue,n)
            end
        end
    end
    outside
end
function boundary_compact(points)
    ps=unique(points)
    length(ps)<3 && return ps
    all(p->abs(cross2(ps[1],ps[2],p))<1e-12,ps) && return ps
    triangles=scattered_triangles(ps)
    radii=[circumradius(ps[a],ps[b],ps[c]) for (a,b,c) in triangles]
    order=sortperm(radii); d=DSU(length(ps)); used=falses(length(ps)); critical=Inf
    nused=0; ncomponents=0
    # MATLAB one-region includes alpha regions touching at a vertex. Requiring
    # a shared triangle edge would enlarge the boundary at such contacts.
    for k in order
        a,b,c=triangles[k]
        for vertex in (a,b,c)
            if !used[vertex]
                used[vertex]=true; nused+=1; ncomponents+=1
            end
        end
        for (u,v) in ((a,b),(b,c))
            if root!(d,u)!=root!(d,v)
                union!(d,u,v); ncomponents-=1
            end
        end
        # Exactly the same all-used, one-component criterion, without scanning
        # every vertex after each triangle on large real event boundaries.
        if nused==length(ps) && ncomponents==1
            critical=radii[k]; break
        end
    end
    # MATLAB boundary's degenerate-spectrum guard selects the convex hull.
    alpha=maximum(radii)-critical<1e-3*maximum(radii) ? Inf : critical
    outside=alpha_exterior(triangles,radii,alpha)
    edgecounts=Dict{Tuple{Int,Int},Int}(); oriented=Dict{Tuple{Int,Int},Tuple{Int,Int}}()
    for (k,(a,b,c)) in enumerate(triangles)
        outside[k] && continue
        for e in ((a,b),(b,c),(c,a))
            key=minmax(e...); edgecounts[key]=get(edgecounts,key,0)+1; oriented[key]=e
        end
    end
    next=Dict{Int,Vector{Int}}()
    for (key,n) in edgecounts
        n==1 || continue; a,b=oriented[key]; push!(get!(next,a,Int[]),b)
    end
    foreach(sort!,values(next))
    loops=Vector{Vector{Tuple{Float64,Float64}}}()
    while !isempty(next)
        # Euler walk retains repeated contact vertices. A single successor per
        # vertex would discard a lobe when two boundary arcs meet at a point.
        start=minimum(keys(next)); stack=Int[start]; circuit=Int[]
        while !isempty(stack)
            v=last(stack)
            if haskey(next,v) && !isempty(next[v])
                w=pop!(next[v]); isempty(next[v]) && delete!(next,v); push!(stack,w)
            else
                push!(circuit,pop!(stack))
            end
        end
        push!(loops,ps[reverse!(circuit)])
    end
    isempty(loops) && return ps
    areas=[abs(sum(p[i][1]*p[mod1(i+1,length(p))][2]-p[i][2]*p[mod1(i+1,length(p))][1] for i in eachindex(p))) for p in loops]
    loops[argmax(areas)] # MATLAB fills all holes before extracting boundary.
end
function scattered_triangles(points)
    connectivity=MiniQhull.delaunay(points)
    triangles=NTuple{3,Int}[]
    sizehint!(triangles,size(connectivity,2))
    for column in eachcol(connectivity)
        a,b,c=Int(column[1]),Int(column[2]),Int(column[3])
        cross2(points[a],points[b],points[c])<0 && ((b,c)=(c,b))
        # Orient counterclockwise for the adaptive containment predicate, then
        # choose a stable cyclic representation for deterministic query order.
        push!(triangles,a<b && a<c ? (a,b,c) : b<c ? (b,c,a) : (c,a,b))
    end
    sort!(triangles)
end
function linear_scattered(points,values,qgridx,qgridy)
    length(points)>=3 || throw(ArgumentError("at least three noncollinear interpolation points required"))
    all(p->abs(cross2(points[1],points[2],p))<1e-12,points) && throw(ArgumentError("collinear interpolation points"))
    out=fill(NaN,size(qgridx))
    triangles=scattered_triangles(points)
    xmin,xmax=extrema(first.(points)); ymin,ymax=extrema(last.(points))
    nbin=clamp(ceil(Int,sqrt(length(triangles)/4)),4,128)
    sx=nbin/(xmax-xmin); sy=nbin/(ymax-ymin)
    bx(x)=clamp(floor(Int,(x-xmin)*sx)+1,1,nbin)
    by(y)=clamp(floor(Int,(y-ymin)*sy)+1,1,nbin)
    bins=[Int[] for _ in 1:nbin^2]
    for (id,t) in enumerate(triangles)
        a,b,c=points[t[1]],points[t[2]],points[t[3]]
        for iy in by(min(a[2],b[2],c[2])):by(max(a[2],b[2],c[2])),
            ix in bx(min(a[1],b[1],c[1])):bx(max(a[1],b[1],c[1]))
            push!(bins[ix+(iy-1)*nbin],id)
        end
    end
    predicates=DT.AdaptiveKernel()
    for p in eachindex(out)
        q=(qgridx[p],qgridy[p])
        xmin<=q[1]<=xmax && ymin<=q[2]<=ymax || continue
        # Candidate indexing avoids the library point locator's repeated global
        # graph scans. Its adaptive containment predicate still decides inside,
        # on-edge, or outside; the Delaunay mesh and interpolation are unchanged.
        triangle_id=0
        for id in bins[bx(q[1])+(by(q[2])-1)*nbin]
            a,b,c=triangles[id]
            if !DT.is_outside(DT.point_position_relative_to_triangle(predicates,points[a],points[b],points[c],q))
                triangle_id=id; break
            end
        end
        triangle_id==0 && continue
        t=triangles[triangle_id]
        a,b,c=t; pa=points[a]; pb=points[b]; pc=points[c]
        den=cross2(pa,pb,pc)
        wa=cross2(q,pb,pc)/den; wb=cross2(pa,q,pc)/den; wc=cross2(pa,pb,q)/den
        # NaNs propagate through a triangle, matching scatteredInterpolant.
        weights=(wa,wb,wc); ids=(a,b,c); value=0.0
        for j in 1:3
            # MATLAB handles exact vertex/edge hits without 0*NaN contamination.
            weights[j]==0 && continue
            value+=weights[j]*values[ids[j]]
        end
        out[p]=value
    end
    out
end
function temporal_normalize(snapshot)
    n=size(snapshot,3); n>=2 || throw(ArgumentError("MATLAB interp1 requires at least two occupied days"))
    out=fill(NaN,100,100,5)
    tb=collect(1:n)./n
    for (j,t) in enumerate(0.0:0.25:1.0)
        lo=clamp(searchsortedlast(tb,t),1,n-1); hi=lo+1; w=(t-tb[lo])/(tb[hi]-tb[lo])
        for k in 1:10000
            a=snapshot[k+(lo-1)*10000]; b=snapshot[k+(hi-1)*10000]
            out[k+(j-1)*10000]=t==tb[lo] ? a : t==tb[hi] ? b : a+w*(b-a)
        end
    end
    out
end

"""
    spn(tracks, lon, lat, anomaly, norm_flag)

Local MATLAB spatial/temporal normalization, returning
`(radius_s, sst_sp, sst_spn_full, x_full, y_full)`. Arithmetic lon/lat centers,
6371 km equirectangular distances without longitude wrapping; each event's
maximum radius; 100×100 polar grid; Delaunay linear interpolation with no
extrapolation; compact one-region alpha boundary (holes filled). `norm_flag==0`
sets values outside the event boundary to zero. Time is interpolated/extrapolated
from (1:n)/n to 0:0.25:1. The final output has size 100×100×5×events×1.

Anomaly starts at the earliest track day, as in the original, and only the first
variable is used for 4-D input. Geometry is separately validated against MATLAB;
cocircular triangulations may select different diagonals for nonlinear data.
"""
function spn(tracks::AbstractVector{Track3D},lon,lat,anomaly::AbstractArray,norm_flag::Real)
    isempty(tracks) && throw(ArgumentError("spn requires at least one event"))
    nx=length(lon); ny=length(lat)
    ndims(anomaly) in (3,4) || throw(DimensionMismatch("anomaly must be lon×lat×time[×variable]"))
    size(anomaly)[1:2]==(nx,ny) || throw(DimensionMismatch("anomaly grid mismatch"))
    data=ndims(anomaly)==4 ? @view(anomaly[:,:,:,1]) : anomaly
    firstday=minimum(minimum(t.day) for t in tracks); lastday=maximum(maximum(t.day) for t in tracks)
    size(data,3)>=lastday-firstday+1 || throw(DimensionMismatch("anomaly must start at earliest event day and cover all events"))
    centers=Vector{Vector{Tuple{Float64,Float64}}}(); radii=Float64[]
    for tr in tracks
        length(tr.day)==length(tr.xloc)==length(tr.yloc) || throw(DimensionMismatch("track fields mismatch"))
        cs=Tuple{Float64,Float64}[]; r=Float64[]
        for j in eachindex(tr.day)
            xs=Int.(tr.xloc[j]); ys=Int.(tr.yloc[j]); isempty(xs) && throw(ArgumentError("empty daily event geometry"))
            c=(nanmean(lon[x] for x in xs),nanmean(lat[y] for y in ys)); push!(cs,c)
            push!(r,maximum(geodist(lat[ys[k]],c[2],lon[xs[k]],c[1])[1] for k in eachindex(xs)))
        end
        push!(centers,cs); push!(radii,maximum(r))
    end
    snapshots=Vector{Array{Float64,3}}(); normalized=fill(NaN,100,100,5,length(tracks),1)
    for (i,tr) in enumerate(tracks)
        x,y=polar_grid(radii[i]); circle=[(x[100,j],y[100,j]) for j in 1:100]
        circle_lookup=PolygonLookup(circle)
        snap=fill(NaN,100,100,length(tr.day))
        for j in eachindex(tr.day)
            cx,cy=centers[i][j]; allpoints=[(geodist(cy,lat[yy],cx,lon[xx])[2:3]) for yy in 1:ny for xx in 1:nx]
            ids=findall(q->inpolygon(q,circle_lookup),allpoints)
            points=allpoints[ids]; values=vec(data[:,:,tr.day[j]-firstday+1])[ids]
            field=linear_scattered(points,values,x,y)
            if norm_flag==0
                pix=[Int(tr.xloc[j][k])+(Int(tr.yloc[j][k])-1)*nx for k in eachindex(tr.xloc[j])]
                border=boundary_compact(allpoints[pix])
                border_lookup=PolygonLookup(border)
                for k in eachindex(field); inpolygon((x[k],y[k]),border_lookup) || (field[k]=0); end
            end
            snap[:,:,j]=field
        end
        push!(snapshots,snap); normalized[:,:,:,i,1]=temporal_normalize(snap)
    end
    xu,yu=polar_grid(1.0)
    radii,snapshots,normalized,xu,yu
end
