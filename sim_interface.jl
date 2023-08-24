# @time includet("landscape_simulation.jl")

tspan = 1630:1:2000
output = MakieOutput(init_state;
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
    axis2 = Axis(fig[2, 1])
    landcover = Observable(Array(frame[].landcover))
    ks = Array(frame[].landcover)
    known_slices = Observable(ks)
    on(frame) do f
        if time[] == 1
            println("zeroing")
            ks .= 0
        end
        landcover[] = f.landcover
        t = tspan[time[]]
        map(NamedTuple(history), NamedTuple(states)) do timeseries, state
            if t in lookup(timeseries, Ti)
                # timeseries = history.native
                # t = 1600
                slice = view(timeseries, Ti(At(t)))
                for I in eachindex(slice)
                    if slice[I]
                        ks[I] = state
                    end
                end
            end
            nothing
        end
        notify(landcover)
        notify(known_slices)
    end
    colormap = cgrad(:Isfahan2, length(states)+1; categorical=true)
    hm = Makie.image!(axis1, landcover; colorrange=(first(states) -1.5, last(states) + 0.5), colormap, interpolate=false)
    Makie.image!(axis2, known_slices; colorrange=(first(states) -1.5, last(states) + 0.5), colormap, interpolate=false)
    ticks = (collect(0:length(states)), vcat(["mask"], collect(string.(propertynames(states)))))
    Colorbar(fig[1, 2], hm; ticks)
    return nothing
end

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
