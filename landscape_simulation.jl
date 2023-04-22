using LandscapeChange
using Distributions
using DimensionalData
using DimensionalData.LookupArrays
using StatsPlots
using Unitful
using StaticArrays
using ReverseStackTraces

includet("raster_common.jl")
includet("roads.jl")
includet("travel_cost.jl")
includet("tabular_data.jl")
includet("/home/raf/.julia/dev/LandscapeChange/src/makie.jl")

println("Getting timeline slices...")
includet("map_file_list.jl")
files = get_map_files()
slices = make_raster_slices(masks, lc_categories)
lc = map(slices.mus.timelines.lc) do A
    reverse(A; dims=Y)
end |> x -> set(x, Ti=>Intervals(End())) |> x -> set(x, Ti=>Irregular((0, 1992)))
# a = @animate for A in lc
#     Plots.plot(A; legend=false, clims=(1, 5))
# end
# Plots.gif(a, "timeseries.gif", fps=1)
# Plots.plot(lc; size=(2000,1300), legend=false, clims=(1, 5))
# savefig("mauritius_timeline.png")

const P = Param

states = NamedVector(native=1, cleared=2, abandoned=3, urban=4, forestry=5)
# to category from category
logic = NamedVector(
    native    = (native=true,  cleared=false, abandoned=false, urban=false, forestry=false,),
    cleared   = (native=true,  cleared=true,  abandoned=true,  urban=false, forestry=false,),
    abandoned = (native=false, cleared=true,  abandoned=true,  urban=false, forestry=false,),
    urban     = (native=false, cleared=false, abandoned=false, urban=true,  forestry=false,),
    forestry  = (native=false, cleared=false, abandoned=false, urban=false, forestry=true,),
)

function can_change(states, logic, to, from, checked=())
    if logic[to][from]
        return true
    else
        return map(states) do s
            if s == to || s in checked
                false
            elseif logic[to][s]
                can_change(states, logic, s, from, (checked..., to))
            else
                false
            end
        end |> any
    end
end
# to category from category with intermediate steps
indirect_logic = map(states) do s1
    map(states) do s2
        can_change(states, logic, s1, s2)
    end
end
pairs(indirect_logic)

function countcats(data, categories=union(data))
    map(categories) do cat
        cat => count(==(cat), data)
    end
end
Plots.plot(human_pop_timeline.reu)
Plots.plot!(human_pop_timeline.mus)

best_slices = lc[At([1600, 1709, 1723, 1772, 1810, 1905, 1992])]
Plots.plot(best_slices)
cat_counts = map(best_slices) do slice
    map(states) do state
        count(x -> x === state, slice)
    end |> NamedTuple
end;
Plots.scatter(human_pop_timeline.mus)
pop_cat_counts = map(lookup(cat_counts, Ti)) do t
    (; cat_counts[At(t)]..., pop=human_pop_timeline.mus[At(t)])
end

using GLM, StatsPlots
# Cleared land is used for urbanisation by 1992, so don't use it in the model
cleared_model = lm(@formula(cleared ~ pop^2 + pop), pop_cat_counts)
urban_model = lm(@formula(urban ~ pop^2 + pop), pop_cat_counts)
ti = dims(human_pop_timeline.mus, Ti)
pops = map(pop -> (; pop), human_pop_timeline.mus)
cleared_pred = DimArray(predict(cleared_model, parent(pops)), ti)
urban_pred = DimArray(predict(urban_model, parent(pops)), ti)
lc_predictions = map((cleared, urban) -> (; cleared, urban), cleared_pred, urban_pred)
Plots.plot(cleared_pred)
Plots.scatter!(map(x -> x.urban, cat_counts))
Plots.plot!(urban_pred)
Plots.scatter!(map(x -> x.cleared, cat_counts))
Plots.plot!(human_pop_timeline.mus)

b = (; bounds=(0.0, 2.0))
inertia = NamedVector(
    native=Param(0.2; b...),
    cleared=Param(0.6; b...),
    abandoned=Param(0.1; b...),
    urban=Param(0.9; b...),
    forestry=Param(0.9; b...),
)
b = (; bounds=(1.0, 2.0))

landscape_events = (
    mus = [
        (year=1638, n=50, geometry=(X=57.7228, Y=-20.3754)), # "First Dutch settlement"
        (year=1664, n=160, geometry=(X=57.7228, Y=-20.3754)),  # "Second dutch settlement"
        (year=1721, n=524, geometry=(X=57.7228, Y=-20.3754)),  # "Colony at grand port founded", "French",
        (year=1721, n=20#=?=#, geometry=(X=57.5012, Y=-20.1597)),  # "Colony at port louie began", "French",
        (year=1735, n=20#=?=#, geometry=(X=57.5012, Y=-20.1597)),  # "Colony at port louie as capital", "French",
    ],
    reu = [
        (year=1665, n=20#=?=#, geometry=(X=57.5012, Y=-20.1597)),  # "Colony at port louie as capital", "French",
    ],
)
fixed = false # or a raster mask of fixed landcover
human_suitability = map(travel_times, slope_stacks, dems) do tt, ss, dem
    slope_suitability = (1 .- ss.slope) .^ 2 # Steep bad, flat good
    slices = map(tt) do travel_time
        travel_suitability = (1 .- travel_time / 30u"hr")
        reverse(replace_missing(travel_suitability .* slope_suitability, 0.0); dims=Y)
    end
    Rasters.combine(slices, Ti)
end

# plot(suitability.mus; clims=(0, 1), legend=false)
# plot(suitability.reu; clims=(0, 1), legend=false)
suitability = map(human_suitability) do hs
    native = ones(axes(hs))
    cleared = hs .^ 2
    abandoned = 1 .- hs
    urban = hs .^ 2
    forestry = ones(axes(hs))
    map(native, cleared, abandoned, urban, forestry) do n, c, a, u, f
        NamedVector(native=n, cleared=c, abandoned=a, urban=u, forestry=f)
    end
end


# The amount of influence neighbors have
b = (; bounds=(1.00, 5.0))
transitions = (
    native = (
        native=Exponential(P(1.7; b..., label="native from native")),
        cleared=0.0,#Exponential(P(1.0; b...)),
        abandoned=0.0,
        urban=0.0,
        forestry=0.0,#Exponential(P(1.0; b...)),
    ),
    cleared = (
        native=0.0,
        cleared=Exponential(P(1.7; b..., label="cleared from cleared")),
        abandoned=0.0,#Exponential(P(1.0; b...)),
        urban=Exponential(P(1.0; b..., label="cleared from urban")),
        forestry=0.0,#Exponential(P(1.0; b...)),
    ),
    abandoned = (
        native=0.0,
        cleared=0.0,#Exponential(P(1.0; b...)),
        abandoned=Exponential(P(1.7; b..., label="abandoned from abandoned")),
        urban=0.0,#Exponential(P(1.0; b...)),
        forestry=0.0,#Exponential(P(1.0; b...)),
    ),
    urban = (
        native=0.0,
        cleared=Exponential(P(1.7; b..., label="urban from cleared")),
        abandoned=0.0,#Exponential(P(1.0; b...)),
        urban=Exponential(P(1.7; b..., label="urban from urban")),
        forestry=0.0,#Exponential(P(1.0; b...)),
    ),
    forestry = (
        native=0.0,
        cleared=0.0,
        abandoned=0.0,
        urban=0.0,
        forestry=Exponential(P(1.7; b..., label="forestry from forestry")),
    ),
)

# The logic of sequential category change - can a transition happen at all
revmasks = map(masks) do A
    reverse(A; dims=Y)
end

# Human Population and species introduction events
eventrule = let events=landscape_events.mus,
                states=states,
                D=dims(revmasks.mus),
                lc=lc
    SetGrid{:landcover,:landcover}() do data, l1, l2
        current_year = currenttime(data)
        if current_year in Rasters.lookup(lc, Ti)
            println("Updating slice...")
            l2 .= lc[At(current_year)]
        end
        for event in events
            if event.year == current_year
                p = event.geometry
                I = DimensionalData.dims2indices(D, (X(Contains(p.X)), Y(Contains(p.Y))))
                Iu = map(i -> i-1:i+1, I)
                Ic = map(i -> i-10:i+10, I)
                l1[Iu...] .= states.urban
                l1[Ic...] .= states.cleared
                l2[Iu...] .= states.urban
                l2[Ic...] .= states.cleared
                # l2 .= states.cleared
            end
        end
    end
end
pressure = NamedVector(
    native=Param(1.0; b...),
    cleared=Param(1.5; b...),
    abandoned=Param(1.0; b...),
    urban=Param(1.0; b...),
    forestry=Param(1.0; b...),
)
pressure = let preds=lc_predictions
    leverage=Param(1.0; bounds=(-50.0, 50.0))
    native=Param(1.0; b...)
    # cleared=Param(1.5; b...),
    abandoned=Param(1.0; b...)
    # urban=Param(1.4; b...)
    forestry=Param(1.0; b...)
    (data, rule) -> begin
        predicted = preds[At(currenttime(data))]
        println(stdout, predicted)
        stategrid = data[:landcover]
        ngridcells = if isnothing(DynamicGrids.mask(data))
            length(stategrid)
        else
            sum(DynamicGrids.mask(data))
        end
        ncleared = sum(==(rule.states.cleared), stategrid)
        nurban = sum(==(rule.states.urban), stategrid)
        println(stdout, (ncleared, nurban))
        urban = leverage * (predicted.urban - nurban) / ngridcells
        cleared = leverage * (predicted.cleared - ncleared) / ngridcells
        NamedVector(; native, cleared, abandoned, urban, forestry)
    end
end

n, c, a, u, f = states
precursors = (
    native =    SA[n, n, n, n],
    cleared =   SA[n, c, a, f],
    abandoned = SA[c, a, a, f],
    urban =     SA[n, c, a, f],
    forestry =  SA[n, c, a, f],
)
staterule = BottomUp{:landcover}(;
    neighborhood=Moore(4),
    states,
    inertia,
    transitions,
    logic=(; precursors, direct=logic, indirect=indirect_logic),
    pressure,
    suitability=Aux{:suitability}(),
    fixed=false,
    perturbation=P(0.1; bounds=(0.0, 1.0), label="perturbation"),
)
ruleset = Ruleset(eventrule, staterule; proc=ThreadedCPU());
landcover = ones(Int, DimensionalData.commondims(suitability.mus, (X, Y)); missingval=0) |>
    A -> mask!(A; with=revmasks.mus)
init_state = (; landcover)
history = Rasters.combine(lc)
# Rasters.rplot(history; colorrange=(1, 5))
aux = (; suitability=suitability.mus, history)
output = ResultOutput(init_state;
    aux,
    tspan=1630:1:1650,
    store=false,
    mask=revmasks.mus,
    boundary=Remove(),
    padval=0,
)
simdata = DynamicGrids.SimData(output, ruleset);
sim!(output, ruleset; simdata, printframe=true);
# set_theme!(backgroundcolor=:white)

# plot(RasterStack(slices.mus.timelines.urban); legend=false, size=(1000, 850), margin=-1mm)
# plot(RasterStack(slices.mus.timelines.lc); legend=false, clims=(0, 5), size=(1000, 850), margin=-1mm)
# plot(last(slices.mus.timelines.lc); legend=false, clims=(0, 5), size=(1000, 850), margin=-1mm)
# savefig("mus_landcover_maps.png")

# lc_slices = collect(slices.mus.timelines.lc)
# plot(last(lc_slices))
# lc_slices
# uncertain_lc = mask(last(lc_slices); with=replace_missing(lc_2017.mus.uncertain, false))
# plot(uncertain_lc)
# countcats(uncertain_lc, lc_categories) |> pairs
# slices.mus.timelines
# lc_years = map(keys(slices.mus.timelines.lc)) do key
#     parse(Int, string(key)[end-3:end])
# end |> collect
# lc_timeline = RasterSeries(lc_slices, Ti(lc_years))
# lc_fractions = map(A -> LandscapeChange.cover_fraction(A; categories=lc_categories), lc_timeline)

# Plots.plot(getproperty.(lc_fractions, :abandoned); labels="abandoned")
# Plots.plot!(getproperty.(lc_fractions, :cleared); labels="cleared")
# Plots.plot!(getproperty.(lc_fractions, :native); labels="native")
# Plots.plot!(getproperty.(lc_fractions, :forestry); labels="forestry")
# Plots.plot!(getproperty.(lc_fractions, :urban); legend=:left, labels="urban", title="Mauritius land cover history")
# Plots.plot!(twinx(), human_pop_timeline.mus; legend=:top, color=:black, labels="human population")
# savefig("mus_landcover_hist.png")

# pop_at_slice = human_pop_timeline.mus[At(lookup(lc_fractions, Ti))]

# table = map(lc_fractions, pop_at_slice) do f, pop
#     NamedTuple(merge(f, (; pop)))
# end |> DataFrame

# f = @formula(x ~ a + b + c + d)
# length(pop_at_slice)
# model = lm(@formula(1 - native ~ pop^2 + pop), table)
# @formula(1 - native ~ pop^2 + pop)
# pops = DataFrame((pop = 1:100:1000000,))
# plot(pops.pop, predict(model, pops))
# Plots.scatter!(table.pop, 1 .- table.native)

# persistence = map(GeometryBasics.TupleView{2,1}(lc_slices)) do xs
#     LandscapeChange.cover_persistence(xs...; categories=lc_categories)
# end
# persistence_timeline =
#     DimArray(persistence, set(dims(lc_timeline, Ti)[begin:end-1], Sampled(; sampling=Intervals(Start()), span=Irregular(first(lc_years), last(lc_years)))))
# annual_persistence = map(enumerate(persistence_timeline)) do (i, t)
#     b = cellbounds(persistence_timeline, i)[1]
#     years = b[2] - b[1]
#     1 .- ((1 .- t) ./ years)
# end
# annual_persistence_timeline = rebuild(persistence_timeline, annual_persistence)

# Plots.plot(getproperty.(annual_persistence_timeline, :abandoned); labels="abandoned")
# Plots.plot!(getproperty.(annual_persistence_timeline, :cleared); labels="cleared")
# Plots.plot!(getproperty.(annual_persistence_timeline, :native); labels="native")
# Plots.plot!(getproperty.(annual_persistence_timeline, :urban); labels="urban", title="Annual persistence of land cover classes", ylims=(0.96, 1.0))
# Plots.plot!(getproperty.(annual_persistence_timeline, :forestry); legend=:left, labels="forestry", title="Annual persistence of land cover classes", ylims=(0.96, 1.0))
# Plots.plot!(twinx(), human_pop_timeline.mus; legend=:bottomleft, color=:black, labels="human population")
# savefig("annual_persistence.png")

# transition_probs = map(GeometryBasics.TupleView{2,1}(lc_slices)) do xs
#     LandscapeChange.cover_change(xs...; categories=lc_categories)
# end
# transition_timeline =
#     DimArray(transition_probs, set(dims(lc_timeline, Ti)[begin:end-1], Sampled(; sampling=Intervals(Start()), span=Irregular(first(lc_years), last(lc_years)))))
# annual_transition = map(enumerate(transition_timeline)) do (i, ts)
#     b = cellbounds(transition_timeline, i)[1]
#     years = b[2] - b[1]
#     map(ts) do t
#         t ./ years
#     end
# end
# annual_transition_timeline = rebuild(transition_timeline, annual_transition)
# annual_transition_timeline[1]
# getproperty.(annual_transition_timeline, :native)
# x = annual_transition_timeline[1]


