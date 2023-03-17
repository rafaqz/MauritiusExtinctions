using LandscapeChange
using Distributions
using DimensionalData
using DimensionalData.LookupArrays
using StatsPlots

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
end
Plots.plot(lc; size=(2000,1300), legend=false, clims=(1, 5))
savefig("mauritius_timeline.png")

const P = Param

function countcats(data, categories=union(data))
    map(categories) do cat
        cat => count(==(cat), data)
    end
end

b = (; bounds=(0.0, 2.0))
inertia = NamedVector(
    native=Param(0.2; b...),
    cleared=Param(0.6; b...),
    abandoned=Param(0.1; b...),
    urban=Param(0.9; b...),
    forestry=Param(0.9; b...),
)
b = (; bounds=(1.0, 2.0))
pressure = NamedVector(
    native=Param(1.0; b...),
    cleared=Param(1.5; b...),
    abandoned=Param(1.0; b...),
    urban=Param(1.4; b...),
    forestry=Param(1.9; b...),
)
states = NamedVector(lc_categories)

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
    cleared = h
    abandoned = 1 .- hs
    urban = hs
    forestry = ones(axes(hs))
    map(native, cleared, abandoned, urban, forestry) do n, c, a, u, f
        NamedVector(native=n, cleared=c, abandoned=a, urban=u, forestry=f)
    end
end

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

revmasks = map(masks) do A
    reverse(A; dims=Y)
end

staterule = BottomUp{:landcover}(;
    neighborhood=Window(5),
    states,
    inertia,
    transitions,
    pressure,
    suitability=Aux{:suitability}(),
    fixed=false,
    perturbation=P(0.1; bounds=(0.0, 1.0), label="perturbation"),
)

lc = map(slices.mus.timelines.lc) do A
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
ruleset = Ruleset(eventrule, staterule; proc=CPUGPU());
landcover = ones(Int, DimensionalData.commondims(suitability.mus, (X, Y)); missingval=0) |>
    A -> mask!(A; with=revmasks.mus)
init_state = (; landcover)
output = ResultOutput(init_state;
    aux=(; suitability=suitability.mus),
    tspan=1600:1:1620,
    store=false,
    mask=revmasks.mus,
    padval=0,
)
simdata = DynamicGrids.SimData(output, ruleset);
sim!(output, ruleset; simdata, printframe=true);

set_theme!(backgroundcolor=:white)
output = MakieOutput(init_state;
    aux=(; suitability=suitability.mus),
    tspan=1637:1:2000,
    store=false,
    mask=revmasks.mus,
    boundary=Remove(),
    padval=0,
    ruleset,
    # sim_kw=(; printframe=true),
) do fig, frame
    axis = Axis(fig[1, 1])
    landcover = Observable(Array(frame[].landcover))
    on(frame) do f
        landcover[] = f.landcover
        notify(landcover)
    end
    colormap = cgrad(:Isfahan2, length(states)+1; categorical=true)
    hm = Makie.heatmap!(axis, landcover; colorrange=(-0.5, 5.5), colormap)
    ticks = (collect(0:length(states)), vcat(["mask"], collect(string.(propertynames(states)))))
    Colorbar(fig[1, 2], hm; ticks)
    # recordframe!(io) # record a new framea
    # record(scene, "test.gif") do io
    #     for i = 1:100
    #         func!(scene)     # animate scene
    #     end
    # end
    return nothing
end
display(output.fig)

using StatsPlots
Plots.plot(Exponential(0.1))

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


includet("/home/raf/.julia/dev/LandscapeChange/src/makie.jl")
using DynamicGrids, GLMakie, Random
# rand!(state, 0:10)

x = "13-26/13-14,17-19/2/M"
x = "13-26/14-19/2/M"
x = "1-3/1,4-5/5/N"
x = "4-7/6-8/10/M"
x = "3,5,7,9,11,15,17,19,21,23-24,26/3,6,8-9,11,14-17,19,24/7/M"
x = "0-3,7-9,11-13,18,21-22,24,26/13,17,20-26/4/M"
x = "1-3/1,4-5/5/N"
x = "4/4/5/M"
x = "5-8/6-7,9,12/4/M"
x = "13-26/13-14,17-19/2/M"
x = "13-26/14-19/2/M"
x = "1,4,8,11,13-26/13-26/5/M"

s, b, i, n = split(x, '/')
neighborhood = n == "M" ? Moore{1,3}() : VonNeumann{1,3}()
start = parse(Int, i)
survive, born = map(split.((s, b), ',')) do vals
    ranges = map(vals) do v
        if occursin("-", v) == 1
            a, b = split(v, '-')
            parse(Int, a) : parse(Int, b)
        else
            parse(Int, v)
        end
    end
    Tuple(vcat(ranges...))
end

function make_init()
    state = zeros(Int, 60, 60, 60)
    # rand!(view(state, 40:50, 40:50, 40:50), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, start))
    # rand!(view(state, 25:30, 25:30, 25:30), (0, 0, 0, 0, 0, 0, 0, 0, start))
    rand!(state, (0, 0, 0, 0, 0, 0, 0, 0, start-3, start - 2, start - 1, start))
    # fill!(state, start)
end

i

state = make_init()
extrainit = Dict(map(i -> Symbol("init", i) => make_init(), 1:10)...)

ca = let survive=survive, born=born, start=start
    Neighbors(neighborhood) do data, hood, state, I
    n = count(>(0), neighbors(hood))
    if state == 0
        n in born ? start : 0
    elseif state == 1
        n in survive ? state : state - 1
    else
        state - 1
    end
end
end
ruleset = Ruleset(ca; proc=CPUGPU())
# output = ResultOutput(init; tspan=1:3)

using GLMakie
using OffsetArrays
A = rand(10, 10)
Makie.heatmap(A)
P = PermutedDimsArray(A, (2, 1))
O = OffsetArray(A, (0:9, 0:9))
Makie.heatmap(O)
using DimensionalData
D = rand(X(10), Y(10))
Makie.heatmap(D)
