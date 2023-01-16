using DynamicGrids
using LandscapeChange
using Distributions

include("tabular_data.jl")

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
lc_categories = NamedVector(;
    native=1, cleared=2, abandoned=3, urban=4,
)

landcover = broadcast(_ -> 0, dems.mus)
counts = broadcast(_ -> zero(lc_categories), landcover)


lc_count_rule = LandCoverCount{:landcover,:counts}(Window(1), counts)
extent = Extent((; landcover, counts); tspan=1:10)
sd = DynamicGrids.SimData(extent, lc_count_rule)
sd = step!(sd)
sd[:counts]
x = parent(sd[:counts])

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

inertia = NamedVector(
    native=1.0,
    cleared=1.0,
    abandonned=1.0,
    urban=1.0,
)

WhiteEngalinWeights(
    transition_potential=
    suitability
    inertia
    active
    fixed
    pertubation
end

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
    soiltypes
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
