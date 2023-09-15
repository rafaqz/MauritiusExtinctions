# @time includet("landscape_simulation.jl")

tspan = 1630:1:2000
output = MakieOutput(init_statae;
    aux, tspan,
    fps=100,
    store=true,
    mask=masks.mus,
    boundary=Remove(),
    padval=0,
    ruleset,
    sim_kw=(; printframe=true),
) do fig, time, frame
    axis1 = Axis(fig[1, 1])
    axis2 = Axis(fig[1, 2])
    landcover = Observable(Array(frame[].landcover))
    known_slices = Observable(view(striped_history, Ti(1)))
    on(frame) do f
        landcover[] = f.landcover
        t = tspan[time[]]::Int
        if hasselection(striped_history, Ti(At(t)))
            known_slices[] = view(striped_history, Ti(At(t)))
            notify(known_slices)
        end
        notify(landcover)
    end
    colormap = cgrad(:batlow, length(states)+1; categorical=true)
    hm = Makie.image!(axis1, landcover; colorrange=(first(states) -1.5, last(states) + 0.5), colormap, interpolate=false)
    Makie.image!(axis2, known_slices; colorrange=(first(states) -1.5, last(states) + 0.5), colormap, interpolate=false)
    ticks = (collect(0:length(states)), vcat(["mask"], collect(string.(propertynames(states)))))
    Colorbar(fig[1, 3], hm; ticks)
    return nothing
end

raster = Raster(CHELSA{BioClim}, :tmax)

# Plots.plot(lc)

# using StatsPlots
# Plots.plot(Exponential(0.1))
# Plots.plot!(Exponential(1.0))
# Plots.plot!(Exponential(2.0))
# Plots.plot!(Exponential(3.0))

# fig = Figure()
# ax = Axis(fig[1, 1])
# Makie.plot!(ax, lc[At(1723)]; colormap=:viridis, colorrange=(0, 5))
# Makie.plot!(ax, lc[At(1772)]; colormap=:viridis, colorrange=(0, 5))
# Makie.plot!(ax, lc[At(1810)]; colormap=:viridis, colorrange=(0, 5))
# Makie.plot!(ax, lc[At(1835)]; colormap=:viridis, colorrange=(0, 5))
# Makie.plot!(ax, lc[At(1854)]; colormap=:viridis, colorrange=(0, 5))
# Makie.plot!(ax, lc[At(1870)]; colormap=:viridis, colorrange=(0, 5))
# Makie.plot!(ax, lc[At(1872)]; colormap=:viridis, colorrange=(0, 5))
# Makie.plot!(ax, lc[At(1905)]; colormap=:viridis, colorrange=(0, 5))
# Makie.plot!(ax, lc[At(1968)]; colormap=:viridis, colorrange=(0, 5))
# Makie.plot!(ax, lc[At(1992)]; colormap=:viridis, colorrange=(0, 5))
# delete!(ax)
