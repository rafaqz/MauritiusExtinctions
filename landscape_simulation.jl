using DynamicGrids
using LandscapeChange

# Human Population and species introduction events
landscape_events = (
    human = (
        # (1638, 50, (-20.3754, 57.7228)), # "First Dutch settlement"
        # (1664, 160, (-20.3754, 57.7228)),  # "Second dutch settlement"
        (1721, (-20.3754, 57.7228), 524),  # "Colony at grand port founded", "French", 
        (1721,  (-20.1597, 57.5012), 20#=?=#),  # "Colony at port louie began", "French", 
        (1735, (-20.1597, 57.5012), ?),  # "Colony at port louie as capital", "French", 
    ),
)


# Land use rules
lc_categories = NamedVector(;
    forrested=1, cleared=2, urban=3,
)

landcover = broadcast(_ -> 0, dems.mus)
counts = broadcast(_ -> zero(lc_categories), landcover)

lc_count_rule = LandCoverCount{:landcover,:counts}(Window(1), counts)
extent = Extent((; landcover, counts); tspan=1:10)
sd = DynamicGrids.SimData(extent, lc_count_rule)
sd = step!(sd)
sd[:counts]
x = parent(sd[:counts])

suitability_parameters = (
    DEM=Aux{:DEM}(),
    slope=Aux{:slope}(),
    to_water=Aux{:to_water}(),
    to_minor_ports=Aux{:to_minor_ports}(),
    to_major_ports=Aux{:to_major_ports}(),
    # to_roads=Aux{:to_roads}(),
)

# interaction_matrix = (
#     forest   = lc -> (forest=Param(0.9),  regrowth=0.0,        cleared=Param(0.1),  settled=0.0)[lc],
#     cleared  = lc -> (forest=Param(0.1),  regrowth=Param(0.1), cleared=Param(0.9),  settled=Param(0.01))[lc],
#     settled  = lc -> (forest=Param(0.01), regrowth=0.0,        cleared=Param(0.01), settled=Param(0.99))[lc],
#     regrowth = lc -> (forest=0.0,         regrowth=Param(0.9), cleared=Param(0.1),  settled=0.0)[lc],
# )

neighbor_matrix = (
    forest   = lc -> (forest=1.0, regrowth=9.0, cleared=0.0, settled=0.0)[lc],
    cleared  = lc -> (forest=Param(0.1),  regrowth=Param(0.2), cleared=Param(0.9),  settled=Param(0.01))[lc],
    settled  = lc -> (forest=Param(0.01), regrowth=0.0,        cleared=Param(0.01), settled=Param(0.99))[lc],
    regrowth = lc -> (forest=0.0,         regrowth=Param(0.9), cleared=Param(0.1),  settled=Param(0.05))[lc],
)

suitability_functions = (;
    forest = p -> 1.0, # Forest is always suitable
    regrowth = p -> 1.0, # Regrowth is always suitable
    cleared = p -> p.to_water * p.slope * p.DEM,
    settled = p -> p.to_minor_ports * p.to_major_ports * p.to_water * p.slope * p.DEM,
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


