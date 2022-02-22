using DynamicGrids,
      Dispersal,
      Rasters,
      Setfield, 
      Shapefile,
      StaticArrays,
      LandscapeChange
      # GeoJSON,
      # Plots,
      # StaticBitArrays

include("data.jl")
include("species_data.jl")


# Land use rules

land_use_category = (;
    forrested=1, cleared=2, urban=3,
)

c = NamedVector(
    forest=1, 
    cleared=2, 
    # settled=3
)

counts = broadcast(_ -> zero(c), landscape)

lu_count_rule = LandCoverCount{:landscape,:counts}(Window(1), c)
extent = Extent((; landscape, counts); tspan=1:10)
sd = DynamicGrids.SimData(extent, lu_count_rule)
sd = step!(sd)
sd[:counts]
x = parent(sd[:counts])

suitability_parameters = (
    riverdist=Aux{:riverdist}(),
    # coastdist=Aux{:distance_to_coast}(),
    elevation=Aux{:elevation}(),
    # steepness=Aux{:steepness}(),
)

# interaction_matrix = (
    # forest=lc -> (forest=0.9, cleared=0.1, settled=0.0)[lc],
    # cleared=lc -> (forest=0.1, cleared=0.9, settled=0.01)[lc],
    # settlement=lc -> (forest=0.0, cleared=0.05, settled=0.95)[lc],
# )

suitability_functions = (;
    forest=p -> 1.0,
    cleared=p -> p.riverdist * p.elevation,
    # settlement=p -> p.riverdist * p.elevation,
)

lu_suitability_rule = LandCoverSuitability{Tuple{},:suitability}(
    functions=suitability_functions, params=()
)

lu_potential_rule = LandCoverPotential{:landscape,:suitability}(;
    categories=c, suitability_parameters, suitability_functions
)


suitability = broadcast(_ -> zero(cf), landscape)
riverdist = rand(dims(suitability))
elevation = rand(dims(suitability))
extent = Extent((; landscape, counts, suitability); 
    aux=(; riverdist, elevation),
    tspan=1:10,
)
# sd = DynamicGrids.SimData(extent, lu_suitability_rule)
# sd = step!(sd)
sd = DynamicGrids.SimData(extent, lu_potential_rule)
sd = step!(sd);
first.(sd[:suitability])


# Species rules

# Stress-based approach
# Pros: generic, can be used with presence/absence models
# Cons: combining stresses is arbitrary
pa_land_use_stress = Cell{Tuple{:N,:LU},:S}() do data, (native_species, lu), I
    land_use_stresses[lu]
end

pa_species_interaction_stress = Cell{Tuple{:N,:S}}() do data, (natives, invasives), I
    map(species, natives) do interactions, native
        sum(map(*, interactions, invasives))
    end
end

pa_human_hunting_stress = Neighbors{:Lu,:S}(Kernel(ExponentialKernel(), Window(5))) do data, hood, lu, I
    ncleared = sum(x for x in hood if x isa categories[:cleared])
    map(*, hunting_susceptibiliy, ncleared)
end

pa_population_loss = Cell{Tuple{:N,:S},:N}() do data, (n, s), I
    n && rand() < s
end

stresses = Combine(pa_land_use_stress, pa_species_interaction_stress, pa_human_hunting_stress)

rules = luc, land_clearing_effect, species_interaction_effect

# Presence/absebce grid initialisation
nspecies = 10
species_presences = (x -> SBitArray(UInt16, false for i in 1:S[1])).suitability

interaction_matrix = SArray(rand(10))
lu_response_matrix = map(_ -> (; forrested=rand(), cleared=rand(), urban=rand()), interaction_matrix)
hunting_susceptibility = SArray(rand(size(interaction_matrix))
landuse_susceptibility = SArray(rand(size(interaction_matrix))
landuse = zeros(Bool, S)


species_pops = (SVector(UInt16, false for i in 1:S[1])


output = ArrayOutput((species=species, land_use=land_use))
sim!(rules)





counters = [0, 0, 0]
countval!(counters, x) = counters[x] += 1  
@btime countval!.(Ref(counters), A)
counters = NamedTuple{map(Symbol, Tuple('a':'j'))}(ntuple(_ -> 0, 10))

@btime countelements($A, $counters)
@code_warntype countelements(A, counters)
@code_native countelements(A, counters)
