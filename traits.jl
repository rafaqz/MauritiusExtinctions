using DataFrames
@time using CSV
using GBIF2
using Chain
using TerminalPager

island_species = @chain CSV.File("/home/raf/PhD/Mauritius/mascarine_species.csv") begin
    DataFrame
    filter(:Binomial=> !ismissing, _)
end
m1 = GBIF2.species_match(island_species.Binomial[1])
taxa = map(island_species.Binomial) do sci
    @show sci
    m = GBIF2.species_match(sci)
    if isnothing(m)
        map(_ -> missing, NamedTuple(m1))
    else
        NamedTuple(m)
    end
end |> DataFrame
island_taxa = hcat(island_species, taxa)
island_taxa.species
names(island_taxa)

occurrence_download(; creator="rafaelschouten@gmail.com", species_Key=2486791)
using GBIF2
occurrence_download("0383754-210914110416597")

@time pan_theria = CSV.read("/home/raf/PhD/Traits/PanTHERIA/ECOL_90_184/PanTHERIA_1-0_WR05_Aug2008.txt", DataFrame; missingstring=["-999.00", "-999"]) 
@time elton_mammals = CSV.read("/home/raf/PhD/Traits/EltonTraits/MamFuncDat.txt", DataFrame)
@time elton_birds = CSV.read("/home/raf/PhD/Traits/EltonTraits/BirdFuncDat.txt", DataFrame)
@time lizards = CSV.File("/home/raf/PhD/Traits/Lizards/Appendix S1 - Lizard data version 1.0.csv") |> DataFrame
@time avonet = CSV.File("/home/raf/PhD/Traits/ELEData/ELEData/TraitData/AVONET1_BirdLife.csv") |> DataFrame
@time combine = CSV.File("/home/raf/PhD/Traits/Combine/COMBINE_archives/trait_data_imputed.csv") |> DataFrame
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
filter(x -> !ismissing(x) && x.Binomial ==(target), island_taxa)
filter(x -> !ismissing(x) && x.Species1 ==(target), avonet)
ilisland_birds = innerjoin(island_taxa, avonet;
    on=:Binomial => :Species1, matchmissing=:notequal, makeunique=true,
)
# island_mammals = innerjoin(island_taxa, elton_mammals; on = :species => :Scientific, matchmissing=:notequal)
# island_birds = innerjoin(island_taxa, elton_birds; on = :species => :Scientific, matchmissing=:notequal)
island_mammals = innerjoin(island_taxa, combine;
    on=:Binomial => :iucn2020_binomial, matchmissing=:notequal, makeunique=true,
)
ilisland_birds = innerjoin(island_taxa, avonet;
    on=:Binomial => :Species1, matchmissing=:notequal, makeunique=true,
)
island_lizards = innerjoin(island_taxa, lizards;
    on=:species => :Binomial, matchmissing=:notequal, makeunique=true
)
# island_lizards = innerjoin(island_taxa, lizards;
    # on=:Binomial => :Binomial, matchmissing=:notequal, makeunique=true
# )

found = vcat(island_lizards.Common, island_birds.Common, island_mammals.Common)
notfound = symdiff(found, island_species.Common)
notfoundtaxa = filter(x -> x.Common in notfound, island_taxa)
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
