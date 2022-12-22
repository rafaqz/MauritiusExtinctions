using DataFrames
@time using CSV
using GBIF2
using Chain
using TerminalPager

mascarene_species = @chain CSV.File("/home/raf/PhD/Mascarines/Tables/mascarine_species.csv") begin
    DataFrame
    filter(:Binomial=> !ismissing, _)
end
m1 = GBIF2.species_match(mascarene_species.Binomial[1])
gbif_species_vec = GBIF2.Species[]

taxa = map(mascarene_species.Binomial) do sci
    println(sci)
    m = GBIF2.species_match(sci)
    if isnothing(m)
        map(_ -> missing, NamedTuple(m1))
    else
        s = GBIF2.species(m)
        push!(gbif_species_vec, s)
        NamedTuple(s)
    end
end
mascarene_taxa = hcat(mascarene_species, taxa)

gbif_species = GBIF2.Table(gbif_species_vec)
CSV.write("gbif_species.csv", gbif_species)

# occurrences = map((mus=:MU, reu=:RE)) do country
#     map(gbif_species) do sp
#         ocs = occurrence_search(sp; country, limit=5000)
#         println(sp.species, " in $country")
#         ocs
#     end
# end;
# occurrence_dfs = map(occurrences) do ocs
#     select!(DataFrame(reduce(vcat, ocs)), Not(:institutionKey))
# end
# map(keys(occurrences), occurrence_dfs) do k, df
#     CSV.write("/home/raf/PhD/Mauritius/Data/Occurrences/$(k)_occurrences.csv", df)
# end
occurrence_dfs = map((mus=:mus, reu=:reu)) do k
    CSV.read("/home/raf/PhD/Mauritius/Data/Occurrences/$(k)_occurrences.csv", DataFrame)
end

x = DataFrames.subset(occurrence_dfs.mus, 
    :species => ByRow(x -> occursin("Psittacula echo", x)))
plot(dems.mus)
scatter!(x.decimalLongitude, x.decimalLatitude)

using Interfaces
@interface ArrayInterface (
    mandatory = (
        size = x -> length(x) == prod(size(x)),
        indices = (
            x -> length(eachindex(IndexLinear(), x)) == length(x),
            x -> length(CartesianIndices(axes(x))) == length(x),
            x -> map((l, c) -> x[l] == x[c], CartesianIndices(axes(x)), eachindex(IndexLinear(), x)) 
        ) 
    ),
    optional = (
        setindex! = x -> false # need to write these...,
        broadcast = x -> false,
    )
)

using Test
function test_array_indices(x)
    lininds = eachindex(IndexLinear(), x)
    carinds = CartesianIndices(axes(x))
    for (li,ci) = zip(lininds,carinds)
        @test x[li] == x[ci]
    end
    length(lininds) == length(carinds) == length(x) == prod(size(x))
end

@interface ArrayInterface (mandatory=(indices=test_array_indices,), optional=())

@implements ArrayInterface Array [zeros(10, 10)]
Interfaces.implements(ArrayInterface, Array)
Interfaces.implements(ArrayInterface{:setindex!}, Array)

iucn_reptiles = CSV.read("/home/raf/PhD/Traits/IUCN data/Reptile IUCN/assessments.csv", DataFrame)
names(iucn_reptiles)

@time pan_theria = CSV.read("/home/raf/PhD/Mauritius/Data/Traits/PanTHERIA/ECOL_90_184/PanTHERIA_1-0_WR05_Aug2008.txt", DataFrame; missingstring=["-999.00", "-999"]) 
@time elton_mammals = CSV.read("/home/raf/PhD/Mauritius/Data/Traits/EltonTraits/MamFuncDat.txt", DataFrame)
@time elton_birds = CSV.read("/home/raf/PhD/Mauritius/Data/Traits/EltonTraits/BirdFuncDat.txt", DataFrame)
@time lizards = CSV.File("/home/raf/PhD/Mauritius/Data/Traits/Lizards/Appendix S1 - Lizard data version 1.0.csv") |> DataFrame
@time avonet = CSV.File("/home/raf/PhD/Mauritius/Data/Traits/ELEData/ELEData/TraitData/AVONET1_BirdLife.csv") |> DataFrame
@time combine = CSV.File("/home/raf/PhD/Mauritius/Data/Traits/Combine/COMBINE_archives/trait_data_imputed.csv") |> DataFrame

x = DataFrames.subset(lizards, :Binomial => ByRow(x -> x === "Leiolopisma mauritiana"))
NamedTuple(pairs(x[1, :]))
names(avonet)
names(combine)
names(pan_theria)
names(elton_birds
names(elton_mammals)
reshape(names(lizards), length(names(lizards)) ÷ 2, 2)
combine.iucn2020_binomial
avonet.Species1
filter(x -> x.Species1 == target, avonet)
target = "Lalage newtoni"
filter(x -> !ismissing(x) && x.Binomial ==(target), mascarene_taxa)
filter(x -> !ismissing(x) && x.Species1 ==(target), avonet)
ilisland_birds = innerjoin(mascarene_taxa, avonet;
    on=:Binomial => :Species1, matchmissing=:notequal, makeunique=true,
)
# island_mammals = innerjoin(mascarene_taxa, elton_mammals; on = :species => :Scientific, matchmissing=:notequal)
# island_birds = innerjoin(mascarene_taxa, elton_birds; on = :species => :Scientific, matchmissing=:notequal)
island_mammals = innerjoin(mascarene_taxa, combine;
    on=:Binomial => :iucn2020_binomial, matchmissing=:notequal, makeunique=true,
)
ilisland_birds = innerjoin(mascarene_taxa, avonet;
    on=:Binomial => :Species1, matchmissing=:notequal, makeunique=true,
)
island_lizards = innerjoin(mascarene_taxa, lizards;
    on=:species => :Binomial, matchmissing=:notequal, makeunique=true
)
# island_lizards = innerjoin(mascarene_taxa, lizards;
    # on=:Binomial => :Binomial, matchmissing=:notequal, makeunique=true
# )

found = vcat(island_lizards.Common, island_birds.Common, island_mammals.Common)
notfound = symdiff(found, island_species.Common)
notfoundtaxa = filter(x -> x.Common in notfound, mascarene_taxa)
pager(sort(notfoundtaxa, :Binomial))


Base.convert(String, ::JSON3.Array) = 

@edit setindex!([Dict{Symbol, String}()], x::JSON3.Object{Vector{UInt8}, SubArray{UInt64, 1, Vector{UInt64}, Tuple{UnitRange{Int64}}, true}}, i1::Int64)

found = subset(elton_birds, "scientific" => ByRow(==(binomial)); skipmissing=true)
if nrow(found) > 0 
    found[1, "BodyMass-Value"]
else
    found = subset(elton_mammals, "Scientific")
    nrow(found) > 0 ? found[1, "BodyMass-Value"] : missing
end
elton_mammals

found = subset(elton_birds, "Scientific" => ByRow(x -> !ismissing(x) && occursin("Zosterops mac", x)); skipmissing=true)

elton_birds

species[!, "weight"] = weight
subset(species, :weight => ByRow(x -> ismissing(x)))
subset(species, :weight => ByRow(x -> !ismissing(x)))
