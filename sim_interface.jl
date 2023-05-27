@time includet("landscape_simulation.jl")
output = MakieOutput(init_state;
    aux,
    fps=100,
    tspan=1630:1:2000,
    store=true,
    mask=revmasks.mus,
    boundary=Remove(),
    padval=0,
    ruleset,
    sim_kw=(; printframe=true),
) do fig, frame
    axis = Axis(fig[1, 1])
    landcover = Observable(Array(frame[].landcover))
    on(frame) do f
        landcover[] = f.landcover
        notify(landcover)
    end
    colormap = cgrad(:Isfahan2, length(states)+1; categorical=true)
    hm = Makie.image!(axis, landcover; colorrange=(-0.5, 5.5), colormap)
    ticks = (collect(0:length(states)), vcat(["mask"], collect(string.(propertynames(states)))))
    Colorbar(fig[1, 2], hm; ticks)
    return nothing
end
Plots.plot(lc)



using StatsPlots
Plots.plot(Exponential(0.1))
Plots.plot!(Exponential(1.0))
Plots.plot!(Exponential(2.0))
Plots.plot!(Exponential(3.0))

fig = Figure()
ax = Axis(fig[1, 1])
Makie.plot!(ax, lc[At(1723)]; colormap=:viridis, colorrange=(0, 5))
Makie.plot!(ax, lc[At(1772)]; colormap=:viridis, colorrange=(0, 5))
Makie.plot!(ax, lc[At(1810)]; colormap=:viridis, colorrange=(0, 5))
Makie.plot!(ax, lc[At(1835)]; colormap=:viridis, colorrange=(0, 5))
Makie.plot!(ax, lc[At(1854)]; colormap=:viridis, colorrange=(0, 5))
Makie.plot!(ax, lc[At(1870)]; colormap=:viridis, colorrange=(0, 5))
Makie.plot!(ax, lc[At(1872)]; colormap=:viridis, colorrange=(0, 5))
Makie.plot!(ax, lc[At(1905)]; colormap=:viridis, colorrange=(0, 5))
Makie.plot!(ax, lc[At(1968)]; colormap=:viridis, colorrange=(0, 5))
Makie.plot!(ax, lc[At(1992)]; colormap=:viridis, colorrange=(0, 5))
delete!(ax)
