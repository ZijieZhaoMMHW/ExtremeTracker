"""Curvilinear polar meshes preserve the actual spn sampling coordinates."""
struct EmbeddedTasmanGIF
    path::String
    title::String
end
function Base.show(io::IO,::MIME"text/html",gif::EmbeddedTasmanGIF)
    print(io,"<figure><figcaption>",gif.title,"</figcaption><img style='max-width:100%' src='data:image/gif;base64,",
        base64encode(read(gif.path)),"'></figure>")
end

function polar_mesh(x,y)
    points=Point2f.(vec(x),vec(y))
    indices=Int[]
    nr,nt=size(x)
    for theta in 1:nt-1, r in 1:nr-1
        a=r+(theta-1)*nr; b=a+1; c=a+nr; d=c+1
        r>1 && append!(indices,(a,b,c)) # omit the degenerate center triangle
        append!(indices,(b,d,c))
    end
    points,indices
end

function polar_panel(parent,geometry,field;title="",limits=(0,5),colormap=:thermal)
    ax=Axis(parent;title,xlabel="x / maximum event radius",ylabel="y / maximum event radius",
        aspect=DataAspect(),backgroundcolor=:lightgray)
    points,indices=geometry
    colors=field isa Observable ? field : vec(field)
    plot=mesh!(ax,points,indices;color=colors,colormap,colorrange=limits,
        nan_color=:lightgray,shading=NoShading,rasterize=1)
    angles=range(0,2pi;length=200)
    lines!(ax,cos.(angles),sin.(angles);color=:gray30,linewidth=1)
    xlims!(ax,-1.03,1.03); ylims!(ax,-1.03,1.03)
    ax,plot
end

function tasman_mask(data,day;transparent=false)
    field=fill(transparent ? Float32(NaN) : 0f0,length(data.lon),length(data.lat))
    for (x,y) in zip(data.track.xloc[day],data.track.yloc[day])
        field[Int(x),Int(y)]=1f0
    end
    field
end

function event_map(parent,data,missing,mask;title="")
    ax=Axis(parent;title,xlabel="Longitude (deg E)",ylabel="Latitude (deg N)",aspect=DataAspect())
    heatmap!(ax,data.lon,data.lat,Float32.(missing);colormap=[:white,:lightgray],colorrange=(0,1))
    heatmap!(ax,data.lon,data.lat,mask;colormap=[:crimson,:crimson],colorrange=(0,1),nan_color=:transparent)
    xlims!(ax,extrema(data.lon)...); ylims!(ax,extrema(data.lat)...)
    ax
end

function record_tasman_gif(data,dates,missing,geometry,result,flag,path)
    day=Observable(1)
    title=lift(j->"Tasman Sea MHW · $(dates[j])",day)
    mask=lift(j->tasman_mask(data,j;transparent=true),day)
    fig=Figure(size=flag===nothing ? (880,430) : (1100,470),fontsize=15)
    event_map(fig[1,1],data,missing,mask;title)
    if flag!==nothing
        color=lift(j->vec(result[2][1][:,:,j]),day)
        title2=lift(j->"norm_flag=$flag · $(dates[j])",day)
        _,plot=polar_panel(fig[1,2],geometry,color;title=title2)
        Colorbar(fig[1,3],plot;label="SST anomaly (°C)",ticks=0:1:5)
        Label(fig[2,:],"Gray: missing input · color scale: 0–5 °C · full-event radius fixed",fontsize=13)
    else
        Label(fig[2,:],"Red: tracked MHW cells · gray: input missing on all 230 days",fontsize=13)
    end
    record(fig,path,1:length(dates);framerate=2,px_per_unit=1) do j
        day[]=j
    end
    path
end
