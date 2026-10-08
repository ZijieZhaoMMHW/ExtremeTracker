"""Save publication/export files and show the actual figure in the notebook."""
function emit_figure(fig, name)
    save(joinpath(FIGURES, name * ".png"), fig; px_per_unit=1.5)
    save(joinpath(FIGURES, name * ".pdf"), fig)
    display(fig)
    nothing
end

const EVENT_COLORS = [:white, :royalblue, :darkorange, :seagreen, :purple,
    :crimson, :goldenrod, :teal, :orchid, :slateblue, :sienna]

function label_panel(parent, data; title="", xlabel="Longitude index", ylabel="Latitude index", maxlabel=maximum(data))
    ax = Axis(parent; title, xlabel, ylabel, aspect=DataAspect(),
        xticks=1:max(1,cld(size(data,1),6)):size(data,1),
        yticks=1:max(1,cld(size(data,2),5)):size(data,2))
    n = max(1, Int(maxlabel))
    palette = [EVENT_COLORS[1]; [EVENT_COLORS[mod1(k, length(EVENT_COLORS)-1)+1] for k in 1:n]]
    heatmap!(ax, 1:size(data,1), 1:size(data,2), data;
        colormap=Makie.Categorical(palette), colorrange=(-0.5,n+0.5))
    ax
end

function strict_overlap_equal(actual, expected)
    length(actual)==length(expected) && all(
        all(getfield(a,k)==getfield(b,k) for k in fieldnames(Track3D))
        for (a,b) in zip(actual,expected))
end

function strict_indexed_equal(actual, expected)
    length(actual)==length(expected) && all(
        a.day==b.day && a.indices==b.indices && a.nvoxels==b.nvoxels
        for (a,b) in zip(actual,expected))
end

function field_error(a,b)
    size(a)==size(b) || error("Different field shapes")
    finite = isfinite.(a) .& isfinite.(b)
    (max_abs=any(finite) ? maximum(abs.(a[finite].-b[finite])) : 0.0,
     nan_mismatch=count(isnan.(a).!=isnan.(b)),
     zero_mismatch=count((a.==0).!=(b.==0)))
end

function read_metrics(path)
    values = Dict{String,String}()
    for line in eachline(path)
        occursin(" = ",line) || continue
        key,value=split(line," = ";limit=2)
        values[key]=value
    end
    values
end
