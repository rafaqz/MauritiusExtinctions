using DynamicGrids,
      GeoJSON,
      Plots,
      Rasters,
      Setfield, 
      Shapefile,
      StaticArrays,
      NearestNeighbors,
      StaticBitArrays

include("data.jl")

land_use_category = (;
    forrested=1,
    cleared=2,
    urban=3,
)

# Grid initialisation
S = 100, 100
species = (SBitArray(UInt16, rand(Bool, 4, 4)) for i in 1:S[1], j in 1:S[2])
interaction_matrix = SArray(rand(10))
lu_response_matrix = map(_ -> (; forrested=rand(), cleared=rand(), urban=rand()), interaction_matrix)
hunting_susceptibility = SArray(rand(size(interaction_matrix))
landuse_susceptibility = SArray(rand(size(interaction_matrix))
landuse = zeros(Bool, S)


land_use_effect = Cell{Tuple{:S,:LU},:S}() do data, (p, l), I
    nh_rings = rings(hood)
    map(interaction_matrix, s) do i, s
        if i == zero(i)
            s
        else
            rand() < i ? false : s
        end
    end
end

species_interaction_effect = Cell{:S}() do data, s, I
    map(interaction_matrix, s) do i, s
        if i == zero(i)
            s
        else
            rand() < i ? false : s
        end
    end
end

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

