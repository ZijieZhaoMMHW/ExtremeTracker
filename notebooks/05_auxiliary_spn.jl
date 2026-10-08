# %% [markdown]
# # 05 · Neighborhood smoothing and event space/time normalization
# `grid_nearby_index`, `convzz`, and `spn` are checked against actual outputs of the original MATLAB implementation, including values and NaN masks.
# Real `spn` shapes come from tracked three-dimensional events; values come from the original mhw_ts field. Its physical meaning has not been independently confirmed, so colorbars use input field value rather than SST anomaly.
# %%
include(isfile("bootstrap.jl") ? "bootstrap.jl" : "notebooks/bootstrap.jl")
# %% [markdown]
# ## `convzz`: neighborhood nanmean at ocean centers
# The center ocean mask determines whether a value is produced; no additional land mask is applied to neighbors. Longitude neighborhoods follow the original 0–360° logic.
# The plots use array indices to avoid misleading rectangular heatmap geometry for uneven longitude coordinates spanning 0/360°.
# %%
aux=matread(joinpath(ROOT,"test","fixtures","auxiliary.mat"))
near=grid_nearby_index(vec(aux["lon_aux"]),vec(aux["lat_aux"]),aux["res"])
smoothed=convzz(aux["data"],near,aux["ocean"])
@test all(near[k]==vecmat(aux["near"][k]) for k in eachindex(near))
@test float_match(smoothed,aux["smooth"])
println("Longitude coordinates: ",vec(aux["lon_aux"]))
fig=Figure(size=(1050,380))
for (c,(name,field)) in enumerate([("Input field",aux["data"]),("Ocean center mask",aux["ocean"]),("convzz nanmean",smoothed)])
    ax=Axis(fig[1,c],title=name,xlabel="Longitude index",ylabel="Latitude index",aspect=DataAspect())
    hm=heatmap!(ax,field;colormap=:viridis,nan_color=:lightgray)
    Colorbar(fig[2,c],hm;vertical=false,label="Value")
end
emit_figure(fig,"05_smoothing")
# %% [markdown]
# ## `spn`: both norm_flag values and all five outputs
# Daily arithmetic coordinate centers; a 6371 km planar distance approximation; one maximum radius for each full event; 100×100 polar sampling.
# `norm_flag=0` sets values outside the event boundary to zero; `norm_flag=1` retains the field. Time nodes `(1:n)/n` are resampled at 0,0.25,0.5,0.75,1, including extrapolation to 0.
# Recompute three real events and check radii, unit coordinates, daily samples, and normalized fields. NaN and zero locations must match exactly.
# %%
references=[matread(joinpath(ROOT,"test","fixtures","spn_real_$flag.mat")) for flag in 0:1]
results=Any[]
spn_errors=NamedTuple[]
@testset "Real event spn MATLAB equivalence" begin
    for (flag,f) in enumerate(references)
        actual=spn(overlap_from_mat(f["tr"]),vec(f["lon"]),vec(f["lat"]),f["anomaly"],flag-1)
        push!(results,actual)
        radius,snapshots,normalized,x,y=actual
        @test isapprox(radius,vecmat(f["radius"]);atol=1e-10,rtol=1e-12)
        @test float_match(x,f["xunit"];atol=2e-15)
        @test float_match(y,f["yunit"];atol=2e-15)
        for i in eachindex(snapshots)
            reference=f["snapshots"][i]
            @test float_match(snapshots[i],reference)
            err=field_error(snapshots[i],reference)
            @test err.nan_mismatch==0 && err.zero_mismatch==0
            push!(spn_errors,(flag=flag-1,event=i,err...))
        end
        @test float_match(normalized,reshape(f["normalized"],size(normalized)))
    end
end
println("Error by event and flag: ",spn_errors)
# %% [markdown]
# ## Daily spatial sampling of one event
# Plot unit x/y sampling points colored by value to show the actual circular sampling geometry. These are polar-grid nodes, not an evenly spaced Cartesian grid; matrix rows/columns are not interpreted as uniform x/y coordinates.
# The gray background indicates missing values that are not plotted. Zero values remain part of the color scale.
# %%
event_id=1
radius,snapshots,normalized,xunit,yunit=results[1]
event_tracks=overlap_from_mat(references[1]["tr"])
days=event_tracks[event_id].day
selected=unique(round.(Int,range(1,length(days);length=3)))
finite_values=vcat([vec(snapshots[event_id][isfinite.(snapshots[event_id])]) for (_,snapshots,_,_,_) in results]...)
limits=extrema(finite_values)
fig=Figure(size=(1150,700))
for flag in 0:1, (c,j) in enumerate(selected)
    field=results[flag+1][2][event_id][:,:,j]; finite=isfinite.(field)
    ax=Axis(fig[flag+1,c],title="norm_flag=$flag · day $(days[j])",xlabel="x / event radius",ylabel="y / event radius",aspect=DataAspect(),backgroundcolor=:lightgray)
    scatter!(ax,xunit[finite],yunit[finite];color=field[finite],colormap=:viridis,colorrange=limits,markersize=3)
    xlims!(ax,-1.03,1.03); ylims!(ax,-1.03,1.03)
end
Colorbar(fig[:,4];colormap=:viridis,limits,label="Input field value")
Label(fig[0,:],@sprintf("Real event %d · maximum radius %.1f km",event_id,radius[event_id]),fontsize=21)
emit_figure(fig,"05_spn_spatial")
# %% [markdown]
# ## Five normalized times and MATLAB residuals
# The top row shows Julia norm_flag=0 values; the bottom shows finite Julia − MATLAB residuals. NaN masks are checked separately and exactly, so missing-value discrepancies cannot be hidden in the residual plot.
# %%
reference=reshape(references[1]["normalized"],size(normalized))
fig=Figure(size=(1300,550))
for j in 1:5
    actual=normalized[:,:,j,event_id,1]; expected=reference[:,:,j,event_id,1]
    finite=isfinite.(actual); both=finite.&isfinite.(expected)
    ax=Axis(fig[1,j],title="Normalized t=$((j-1)/4)",xlabel="x / radius",ylabel="y / radius",aspect=DataAspect(),backgroundcolor=:lightgray)
    scatter!(ax,xunit[finite],yunit[finite];color=actual[finite],colormap=:viridis,colorrange=limits,markersize=2.5)
    ax2=Axis(fig[2,j],title="Julia − MATLAB",xlabel="x / radius",ylabel="y / radius",aspect=DataAspect(),backgroundcolor=:lightgray)
    scatter!(ax2,xunit[both],yunit[both];color=actual[both].-expected[both],colormap=:balance,colorrange=(-1e-14,1e-14),markersize=2.5)
    xlims!(ax,-1.03,1.03); ylims!(ax,-1.03,1.03)
    xlims!(ax2,-1.03,1.03); ylims!(ax2,-1.03,1.03)
end
Colorbar(fig[1,6];colormap=:viridis,limits,label="Input field value")
Colorbar(fig[2,6];colormap=:balance,limits=(-1e-14,1e-14),label="Residual")
emit_figure(fig,"05_spn_normalized_and_residual")
