using DynamicGrids,
      Dispersal,
      # GeoJSON,
      # Plots,
      Rasters,
      Setfield, 
      Shapefile,
      StaticArrays
      # StaticBitArrays

include("data.jl")

land_use_category = (;
    forrested=1, cleared=2, urban=3,
)

# Grid initialisation
S = 100, 100
species = (SBitArray(UInt16, rand(Bool, 4, 4)) for i in 1:S[1], j in 1:S[2])
interaction_matrix = SArray(rand(10))
lu_response_matrix = map(_ -> (; forrested=rand(), cleared=rand(), urban=rand()), interaction_matrix)
hunting_susceptibility = SArray(rand(size(interaction_matrix))
landuse_susceptibility = SArray(rand(size(interaction_matrix))
landuse = zeros(Bool, S)


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
output = ArrayOutput((species=species, land_use=land_use))
sim!(rules)

counters = [0, 0, 0]
countval!(counters, x) = counters[x] += 1  
@btime countval!.(Ref(counters), A)
counters = NamedTuple{map(Symbol, Tuple('a':'j'))}(ntuple(_ -> 0, 10))

@btime countelements($A, $counters)
@code_warntype countelements(A, counters)
@code_native countelements(A, counters)
