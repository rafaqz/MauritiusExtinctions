using LandscapeChange
using Distributions
using GeometryBasics
using DimensionalData
using DimensionalData.LookupArrays
using GLM

# using Tyler, TileProviders
# using GLMakie
# tyler = Tyler.Map(Rect2f(57.0, -20.7, 1, 1); provider=Google(), resolution=(1500, 1500))

includet("common.jl")
includet("tabular_data.jl")
includet("raster_common.jl")
includet("map_file_list.jl")

files = get_map_files()
masks = map(boolmask, dems)
lc_categories = NamedVector(native=1, forestry=2, cleared=3, abandoned=4, urban=5)
slices = make_raster_slices(masks, lc_categories)
# plot(RasterStack(slices.mus.timelines.urban); legend=false, size=(1000, 850), margin=-1mm)
# plot(RasterStack(slices.mus.timelines.lc); legend=false, clims=(0, 5), size=(1000, 850), margin=-1mm)
# plot(last(slices.mus.timelines.lc); legend=false, clims=(0, 5), size=(1000, 850), margin=-1mm)
# savefig("mus_landcover_maps.png")

function countcats(data, categories=union(data))
    map(categories) do cat
        cat => count(==(cat), data)
    end
end

plot(last(lc_slices))
uncertain_lc = mask(last(lc_slices); with=replace_missing(lc_2017.mus.uncertain, false))
plot(uncertain_lc)
countcats(uncertain_lc, lc_categories) |> pairs

lc_slices = collect(slices.mus.timelines.lc)
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

# Human Population and species introduction events
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

# Land use rules

landcover = broadcast(_ -> 0, dems.mus)
counts = broadcast(_ -> zero(lc_categories), landcover)

inertia = NamedVector(
    native=1.0,
    cleared=1.0,
    abandoned=1.0,
    urban=1.0,
)
states = NamedVector(
    native=0x01,
    cleared=0x02,
    abandoned=0x03,
    urban=0x03,
)

fixed = false # or a raster mask of fixed landcover
suitability = 1.0 # or a suitability raster or a dynamic Grid?
transition_raster = map(_ -> LandscapeChange.weu_zero(), dems.mus)
state_raster = map(_ -> zero(UInt8), dems.mus)
transition_possibility = (
    forest =     (forest=true,  cleared=true,  abandoned=false, settled=true),
    abandoned = (forest=false, cleared=true,  abandoned=true,  settled=true),
    cleared =    (forest=false, cleared=true,  abandoned=true,  settled=true),
    settled =    (forest=false, cleared=false, abandoned=true,  settled=true),
)

ncleared = 
nforested = 
nabandoned = 
nurban = 
cell_counts = map(ncleared, nforested, nabandoned, nurban) do c, f, a, u
    (cleared=c, forested=f, abandoned=a, urban=u)
end
targets = 

const P = Param
transition_potential = (
    forest = (
        forest=(0.9, 0.0),
        cleared=P(0.1),
        abandoned=0.0,
        settled=0.0,
    ),
    abandoned = (
        forest=(-100.0, 0.0),
        cleared=(P(-0.1), P(0.0)),
        abandoned=(P(0.9), P(0.0)),
        settled=(P(0.05), P(0.1)),
    ),
    cleared = (
        forest=(-1.0, 0.0),
        cleared=(P(0.9), P(0.2)),
        abandoned=(P(0.2), P(0.0)),
        settled=(P(80.0), P(0.0)),
    ),
    settled = (
        forest=(-1.0, 0.0),
        cleared=(P(0.9), P(0.0)),
        abandoned=(P(0.99), P(0.0)),
        settled=(P(0.99), P(-0.2)),
    ),
)

stripparams(transition_potential)

fill_steps(

weight_rule = WhiteEngalinUljeeWeights{:state,:weights}(
    inertia,
    transition_potential=
    suitability=Aux{:suitability}()
    fixed=false,
    pertubation=0.1,
)
update_rule = WhiteEngalinUljeeUpdate{:weights,:state}(Aux(:targets))


# lc_count_rule = LandCoverCount{:landcover,:counts}(Window(1), counts)
# extent = Extent((; landcover, counts); tspan=1:10)
# sd = DynamicGrids.SimData(extent, lc_count_rule)
# sd = step!(sd)
# sd[:counts]
# x = parent(sd[:counts])

exponential1(x) = pdf(Exponential(1), x)
chi5(x) = pdf(Chisq(4), x)

choice_stacks = map(dems, slope_stacks, travel_times) do dem, slope, travel
    RasterStack((
        elevation=dem,
        slope=slope[:slope],
        travel,
    ))
end
plot(choice_stacks.reu)
plot(choice_stacks.mus; c=reverse(colors1))
plot(choice_stacks.mus; c=(colors1))

choice_parameters = (
    elevation=(f=exponential1, scalar=Param(0.001)),
    slope=(f=exponential1, scalar=Param(4.0)),
    travel=(f=exponential1, scalar=Param(0.03)),
    # to_roads=(f=exponential1, scalar=Param(0.03)),
    # to_water=(f=exponential1, scalar=Param(0.01)),
    # to_coast=(f=exponential1, scalar=Param(0.01)),
    # to_minor_ports=(f=exponential1, scalar=Param(0.01)),
    # to_major_ports=(f=exponential1, scalar=Param(0.001)),
)


choice_stacks_scaled = map(choice_stacks) do stack
    map(stack, choice_parameters) do rast, p
        broadcast(rast) do x
            ismissing(x) ? missing : p.f(x * p.scalar)
        end
    end |> RasterStack
end

plot(choice_stacks_scaled.mus; c=:viridis)
plot(choice_stacks_scaled.reu; c=:viridis)
plot(choice_stacks_scaled.mus[:to_water]; c=:viridis)

c = choice_stacks_scaled.mus |> NamedTuple
plot(c.to_roads .* c.to_water .* c.to_major_ports .* c.to_minor_ports .* c.slope; c=:viridis)

choice_aux = map(Aux, namedkeys(choice_parameters))

# interaction_matrix = (
#     forest   = lc -> (forest=Param(0.9),  regrowth=0.0,        cleared=Param(0.1),  settled=0.0)[lc],
#     cleared  = lc -> (forest=Param(0.1),  regrowth=Param(0.1), cleared=Param(0.9),  settled=Param(0.01))[lc],
#     settled  = lc -> (forest=Param(0.01), regrowth=0.0,        cleared=Param(0.01), settled=Param(0.99))[lc],
#     regrowth = lc -> (forest=0.0,         regrowth=Param(0.9), cleared=Param(0.1),  settled=0.0)[lc],
# )

neighbor_matrix = (
    forest   = (forest=1.0,        regrowth=0.0,        cleared=0.0,         settled=0.0),
    cleared  = (forest=Param(0.1), regrowth=Param(0.2), cleared=Param(0.9),  settled=Param(0.01)),
    settled  = (forest=0.0,        regrowth=0.0,        cleared=Param(0.01), settled=Param(0.99)),
    regrowth = (forest=0.0,        regrowth=Param(0.9), cleared=Param(0.1),  settled=Param(0.05)),
)

suitability_functions = (;
    forest = p -> 1.0, # Forest is always suitable
    regrowth = p -> 1.0, # Regrowth is always suitable
    cleared = p -> p.to_water * p.slope * p.DEM * p.soil[:cleared],
    settled = p -> p.travel * p.to_water * p.slope * p.DEM,
)

lu_suitability_rule = LandCoverSuitability{Tuple{},:suitability}(
    functions=suitability_functions, params=suitability_parameters
)

lu_potential_rule = LandCoverPotential{:landcover,:suitability}(;
    categories=lc_categories, suitability_parameters, suitability_functions
)

aux = NamedTuple(
    dems.mus,
    sloperasters.mus,
    distance_to_water.mus,
    distance_to_ports.mus.minor,
    distance_to_ports.mus.major,
    soiltypes,
)

extent = Extent((; landcover, counts, suitability);
    aux, tspan=1600:2000,
)

# sd = DynamicGrids.SimData(extent, lu_suitability_rule)
# sd = step!(sd)
sd = DynamicGrids.SimData(extent, lu_potential_rule)
sd = step!(sd);
first.(sd[:suitability])

using StatsPlots
using Distributions
plot(Exponential(1))
pdf(Exponential(1), 10 * 0.1)

mode = RandomConstraintMatch()
sim(mode, init, target, categories)
