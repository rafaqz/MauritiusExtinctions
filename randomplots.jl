

using ProfileView
using Cthulhu
@profview for i in 1:100000 DynamicGrids.descendable(simdata) end
@descend DynamicGrids.descendable(simdata)
@profview sim!(result_output, ruleset; simdata, proc=CPUGPU()); 

plot(RasterStack(slices.mus.timelines.urban); legend=false, size=(1000, 850), margin=-1mm)
plot(RasterStack(slices.mus.timelines.lc); legend=false, clims=(0, 5), size=(1000, 850), margin=-1mm)
plot(last(slices.mus.timelines.lc); legend=false, clims=(0, 5), size=(1000, 850), margin=-1mm)
savefig("mus_landcover_maps.png")

lc_slices = collect(slices.mus.timelines.lc)
plot(last(lc_slices))
lc_slices
uncertain_lc = mask(last(lc_slices); with=replace_missing(lc_2017.mus.uncertain, false))
plot(uncertain_lc)
countcats(uncertain_lc, lc_categories) |> pairs
slices.mus.timelines
lc_years = map(keys(slices.mus.timelines.lc)) do key
    parse(Int, string(key)[end-3:end])
end |> collect
lc_timeline = RasterSeries(lc_slices, Ti(lc_years))
lc_fractions = map(A -> LandscapeChange.cover_fraction(A; categories=lc_categories), lc_timeline)

Plots.plot(getproperty.(lc_fractions, :abandoned); labels="abandoned")
Plots.plot!(getproperty.(lc_fractions, :cleared); labels="cleared")
Plots.plot!(getproperty.(lc_fractions, :native); labels="native")
Plots.plot!(getproperty.(lc_fractions, :forestry); labels="forestry")
Plots.plot!(getproperty.(lc_fractions, :urban); legend=:left, labels="urban", title="Mauritius land cover history")
Plots.plot!(twinx(), human_pop_timeline.mus; legend=:top, color=:black, labels="human population")
savefig("mus_landcover_hist.png")

pop_at_slice = human_pop_timeline.mus[At(lookup(lc_fractions, Ti))]

table = map(lc_fractions, pop_at_slice) do f, pop
    NamedTuple(merge(f, (; pop)))
end |> DataFrame

f = @formula(x ~ a + b + c + d)
length(pop_at_slice)
model = lm(@formula(1 - native ~ pop^2 + pop), table)
@formula(1 - native ~ pop^2 + pop)
pops = DataFrame((pop = 1:100:1000000,))
plot(pops.pop, predict(model, pops))
Plots.scatter!(table.pop, 1 .- table.native)

persistence = map(GeometryBasics.TupleView{2,1}(lc_slices)) do xs
    LandscapeChange.cover_persistence(xs...; categories=lc_categories)
end
persistence_timeline =
    DimArray(persistence, set(dims(lc_timeline, Ti)[begin:end-1], Sampled(; sampling=Intervals(Start()), span=Irregular(first(lc_years), last(lc_years)))))
annual_persistence = map(enumerate(persistence_timeline)) do (i, t)
    b = cellbounds(persistence_timeline, i)[1]
    years = b[2] - b[1]
    1 .- ((1 .- t) ./ years)
end
annual_persistence_timeline = rebuild(persistence_timeline, annual_persistence)

Plots.plot(getproperty.(annual_persistence_timeline, :abandoned); labels="abandoned")
Plots.plot!(getproperty.(annual_persistence_timeline, :cleared); labels="cleared")
Plots.plot!(getproperty.(annual_persistence_timeline, :native); labels="native")
Plots.plot!(getproperty.(annual_persistence_timeline, :urban); labels="urban", title="Annual persistence of land cover classes", ylims=(0.96, 1.0))
Plots.plot!(getproperty.(annual_persistence_timeline, :forestry); legend=:left, labels="forestry", title="Annual persistence of land cover classes", ylims=(0.96, 1.0))
Plots.plot!(twinx(), human_pop_timeline.mus; legend=:bottomleft, color=:black, labels="human population")
savefig("annual_persistence.png")

transition_probs = map(GeometryBasics.TupleView{2,1}(lc_slices)) do xs
    LandscapeChange.cover_change(xs...; categories=lc_categories)
end
transition_timeline =
    DimArray(transition_probs, set(dims(lc_timeline, Ti)[begin:end-1], Sampled(; sampling=Intervals(Start()), span=Irregular(first(lc_years), last(lc_years)))))
annual_transition = map(enumerate(transition_timeline)) do (i, ts)
    b = cellbounds(transition_timeline, i)[1]
    years = b[2] - b[1]
    map(ts) do t
        t ./ years
    end
end
annual_transition_timeline = rebuild(transition_timeline, annual_transition)
annual_transition_timeline[1]
getproperty.(annual_transition_timeline, :native)
x = annual_transition_timeline[1]


Rasters.rplot(history.abandoned)
Rasters.rplot(history.cleared)
Rasters.rplot(history.urban)
Rasters.rplot(history.native)
Rasters.rplot(history.forestry)

Rasters.rplot(slices.mus.files.landcover_1965.grouped.cleared)


 set_theme!(theme_dark())

 ps = map(lookup(lc, Ti)) do t
     r = replace_missing(lc[Ti(At(t))], NaN)
     figure = Figure(; 
         backgroundcolor=:transparent,
     )
     axis = Axis(figure[1, 1];
         xticklabelsvisible=false, 
         yticklabelsvisible=false,
         xticksvisible=false, 
         yticksvisible=false,
         xgridvisible=false,
         ygridvisible=false,
         bottomspinevisible=false,
         topspinevisible=false,
         leftspinevisible=false,
         rightspinevisible=false,
         aspect=DataAspect(),
     )
     text = "$t"
     textpos = 57.3, -20.0
     Makie.text!(axis, textpos...; text, fontsize=80, align=(:left, :top)) 
     p = Makie.heatmap!(axis, r; 
         colormap=:batlow, 
         colorrange=(0, 6),
     )jj\zz
     save("images/mus_lc_$t.png", figure)
     p
 end

 CairoMakie.activate!()
 GLMakie.activate!()
 savefig("mauritius_timeline.png")
 Rasters.rplot(Rasters.combine(slices.mus.timelines.lc, Ti); 
     colormap=:batlow, colorrange=(0, 5),
     aspect=DataAspect(),
     axis = (
         xticklabelsvisible=false, 
         yticklabelsvisible=false,
         xticksvisible=false, 
         yticksvisible=false,
         xgridvisible=false,
         ygridvisible=false,
         bottomspinevisible=false,
         topspinevisible=false,
         leftspinevisible=false,
         rightspinevisible=false,
     ),
 )
 lc = map(slices.mus.timelines.lc) do A
     reverse(A; dims=Y)
 end |> x -> set(x, Ti=>Intervals(End())) |> x -> set(x, Ti=>Irregular((0, 1992)))
 a = @animate for A in lc
     Plots.plot(A; legend=false, clims=(1, 5))
 end
 Plots.gif(a, "timeseries.gif", fps=1)
 Plots.plot(lc; size=(2000,1300), legend=false, clims=(1, 5))
 savefig("mauritius_timeline.png")

