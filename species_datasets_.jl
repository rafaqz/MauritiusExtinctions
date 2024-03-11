using DynamicGrids, Dispersal, LandscapeChange, ModelParameters
using StaticArrays 
using Statistics
using StatsBase
using Shapefile
using CSV, DataFrames, XLSX, TerminalPager
using GBIF2
includet("raster_common.jl")

# IUCNRedList.set_token("486f559a61285ba396234fc186897b94eda1bd15aa4216a8e1d9f4a8cf40d4c7")
# spec = species_by_category("EX") 
# extinct_species = map(spec["result"]) do row
#     (
#         scientific_name = row["scientific_name"],
#         subspecies      = row["subspecies"],
#         rank            = row["rank"],
#         subpopulation   = row["subpopulation"],
#         taxonid         = row["taxonid"],
#     )
# end
# spec = species_by_category("EW") 
# DataFrame(spec)
# extinct_wild = map(spec["result"]) do row
#     (
#         scientific_name = row["scientific_name"],
#         subspecies      = row["subspecies"],
#         rank            = row["rank"],
#         subpopulation   = row["subpopulation"],
#         taxonid         = row["taxonid"],
#     )
# end

# narratives = map(extinct_species) do s
#     @show s.taxonid
#     species_narrative(s.taxonid)["result"][1]
# end
# narratives_wild = map(extinct_wild) do s
#     @show s.taxonid
#     species_narrative(s.taxonid)["result"][1]
# end
# function nar_rows(narratives) 
#     map(narratives) do row 
#         (
#           habitat              = get(row, "habitat", missing),
#           threats              = get(row, "threats", missing),
#           population           = get(row, "population", missing),
#           conservationmeasures = get(row, "conservationmeasures", missing),
#           rationale            = get(row, "rationale", missing),
#           geographicrange      = get(row, "geographicrange", missing),
#           populationtrend      = get(row, "populationtrend", missing),
#           usetrade             = get(row, "usetrade", missing),
#           taxonomicnotes       = get(row, "taxonomicnotes", missing),
#           species_id           = get(row, "species_id", missing),
#         )
#     end |> DataFrame
# end
# narratives_wild_rows = nar_rows(narratives_wild)
# extinct_wild_all = leftjoin(DataFrame(extinct_wild), narratives_wild_rows; on=:taxonid=>:species_id)

# extinct_all = leftjoin(extinct_df, narrative_rows; on=:taxonid=>:species_id)
# extinct_all.status .= "EX"
# extinct_wild_all.status .= "EW"
# ex_ew = vcat(extinct_all, extinct_wild_all)
# map(enumerate(eachrow(ex_ew)[1:2])) do (i, row)
#     @show i row.scientific_name
#     species_synonyms(row.scientific_name)
# end
# CSV.write("/home/raf/Data/Extinction/redlist/extinct_species.csv", ex_ew; transform=(col, val) -> something(val, missing))

function set_gbif_species!(df, specieskey)
    if !("GBIFSpecies" in names(df))
        df.GBIFSpecies .= ""
    end
    specvec = getproperty(df, specieskey)
    for i in eachindex(specvec) 
        df.GBIFSpecies[i] == "" || continue
        sp = specvec[i]
        ismissing(sp) && continue
        @show sp
        match = GBIF2.species_match(sp)
        isnothing(match) || ismissing(match.species) && continue
        df.GBIFSpecies[i] = match.species 
    end
end

mascarene_species_csv = "/home/raf/PhD/Mascarenes/MauritiusExtinctions/tables/mascarene_species.csv"
@async run(`libreoffice $mascarene_species_csv`)
s = CSV.read(mascarene_species_csv, DataFrame)# |> 
    # x -> subset(x, :Species => ByRow(!ismissing), :GBIFSpecies => ByRow(!ismissing))

# s.GBIFSpecies = copy(s.Species)
# for i in eachindex(s.Species) 
#     sp = s.Species[i]
#     ismissing(sp) && continue
#     match = species_match(sp)
#     isnothing(match) && continue
#     s.GBIFSpecies[i] = match.species 
# end
# s.GBIFSpecies

# CSV.write(mascarene_species_csv, s)
iucn_extinct = CSV.read("/home/raf/Data/Extinction/redlist/extinct_species.csv", DataFrame)
# iucn_extinct.threats |> pager
iucn_csv = "/home/raf/Data/Extinction/redlist/redlist_species_data_39da78ce-d594-4968-8043-489f2765d687/assessments_gbif.csv"
iucn_df = CSV.read(iucn_csv, DataFrame)
set_gbif_species!(iucn_df, :scientificName)
CSV.write(iucn_csv, iucn_df)

# run(`libreoffice /home/raf/Data/Extinction/redlist/redlist_species_data_39da78ce-d594-4968-8043-489f2765d687/assessments_gbif.csv`)

# s_iucn = leftjoin(s, iucn_df; on=:GBIFSpecies, matchmissing=:notequal, makeunique=true)
# s_iucn
# .habitat
# sort!(s_iucn, :GBIFSpecies)
# sort!(s, :GBIFSpecies)
# s.iucn_habitat = s_iucn.habitat 
# s.iucn_rationale = s_iucn.rationale 
# s.iucn_threats = s_iucn.threats 
# joined[!, [:Species, :threats]] |> pager
# broadcast(string, joined.Species, Ref(": "), joined.threats) |> pager


# sp_nothing = filter(x -> isnothing(x[2]), sp_pairs)
# sp = filter(x -> !isnothing(x[2]), sp_pairs)
# gbif_sp = DataFrame(last.(sp)
# first.(sp) .=> gbif_sp.species
# first.(sp)[first.(sp) .!== gbif_sp.species]
# redlist_extinct_csv = "/home/raf/Data/Traits/redlist_species_data_1f74a1f8-0b29-4567-9766-046807e966ca/taxonomy.csv"
# redlist = CSV.read(redlist_extinct_csv, DataFrame; normalizenames=true)

# pantheria_csv = "/home/raf/Data/Traits/PanTHERIA/ECOL_90_184/PanTHERIA_1-0_WR05_Aug2008.txt"
# pantheria = CSV.read(pantheria_csv, DataFrame; normalizenames=true, quoted=false)
# set_gbif_species!(pantheria, :MSW05_Binomial)
pantheria_csv = "/home/raf/Data/Traits/PanTHERIA/ECOL_90_184/PanTHERIA_1-0_WR05_Aug2008_gbif.csv"
CSV.write(pantheria_csv, pantheria)
pantheria_mass = pantheria[!, [:GBIFSpecies, :AdultBodyMass_g]]
s_pantheria = leftjoin(s, pantheria_mass; on=:GBIFSpecies, matchmissing=:notequal, makeunique=true)
sort!(s_pantheria, :GBIFSpecies)

avonet_csv = "/home/raf/Data/Traits/Avonet/ELEData/ELEData/TraitData/AVONET1_BirdLife_gbif.csv"
avonet = CSV.read(avonet_csv, DataFrame; normalizenames=true)
filter(:GBIFSpecies => s -> !ismissing(s) && contains(s, "Crypt"), avonet)
# set_gbif_species!(avonet, :Species1)
# CSV.write(avonet_csv, avonet)
avonet_mass = avonet[!, [:GBIFSpecies, :Mass]]
s_avonet = unique(leftjoin(s, avonet; on=:GBIFSpecies, matchmissing=:notequal, makeunique=true), 3)

filter(r -> r.Species == "Gallinula chloropus", s_avonet)

lizzard_csv = "/home/raf/Data/Traits/Lizards/Appendix S1 - Lizard data version 1.0.csv"
lizzard = CSV.read(lizzard_csv, DataFrame; normalizenames=true)

mass_cols = ["Binomial", "mass_equation_Feldman_et_al_2016_unless_stated_", "intercept", "slope"]
lizzard_mass = lizzard[!, mass_cols]
filter(:Binomial => g -> !ismissing(g) && g == "Proscelotes arnoldi", lizzard_mass)
s_lizzard = leftjoin(s, lizzard_mass; on=:GBIFSpecies, matchmissing=:notequal, makeunique=true)
lizzard_end = filter(r -> !ismissing(r.intercept) && r.Origin == "Endemic", s_lizzard)
# lizzard_end |> pager

elton_bird_csv = "/home/raf/Data/Traits/EltonTraits/BirdFuncDat_gbif.txt"
elton_bird = CSV.read(elton_bird_csv, DataFrame; normalizenames=true)
set_gbif_species!(elton_bird, :Scientific)
CSV.write(elton_bird_csv, elton_bird)
names(elton_bird)

# elton_mammal_csv = "/home/raf/Data/Traits/EltonTraits/MamFuncDat_csv.txt"
# set_gbif_species!(elton_mammal, :Scientific)
elton_mammal_csv = "/home/raf/Data/Traits/EltonTraits/MamFuncDat_gbif.csv"
elton_mammal = CSV.read(elton_mammal_csv, DataFrame; normalizenames=true)
# CSV.write(elton_mammal_csv, elton_mammal)
names(elton_mammal)
elton_mass = vcat(elton_mammal[!, [:GBIFSpecies, :BodyMass_Value]], elton_bird[!, [:GBIFSpecies, :BodyMass_Value]])
s_elton = leftjoin(s, elton_mass; on=:GBIFSpecies, matchmissing=:notequal)

# bird_mass_path = "/home/raf/PhD/Mascarenes/Tables/Bird Mass filled (Jan 22 2015)_WDK.xlsx"
# bird_mass_xl = XLSX.readxlsx(bird_mass_path)
# bird_mass_df = DataFrame(XLSX.eachtablerow(bird_mass_xl["Bird Mass filled"]))
# set_gbif_species!(bird_mass_df, :BirdLife_SpecName)
# CSV.write(bird_mass_csv, bird_mass_df)
bird_mass_csv = "/home/raf/PhD/Mascarenes/Tables/Bird Mass filled (Jan 22 2015)_WDK_gbif.csv"
bird_mass_df = CSV.read(bird_mass_csv, DataFrame)
s_bird_mass = leftjoin(s, bird_mass_df; on=:GBIFSpecies, matchmissing=:notequal, makeunique=true, order=:left)
mean(skipmissing(filter(r -> r.Genus == "Anas", bird_mass_df).MeanBodyMass))
filter(:GBIFSpecies => s -> s === "Cryptopsar", bird_mass_df)


# reptile_mass_path = "/home/raf/PhD/Mascarenes/Tables/Reptile body mass database Meiri 2010.xlsx"
# reptile_mass_xl = XLSX.readxlsx(reptile_mass_path)
# reptile_mass_df = DataFrame(XLSX.eachtablerow(reptile_mass_xl[1]))
# set_gbif_species!(reptile_mass_df, :Name)
# CSV.write(reptile_mass_csv, reptile_mass_df)
reptile_mass_csv = "/home/raf/PhD/Mascarenes/Tables/Reptile body mass database Meiri 2010_gbif.csv"
reptile_mass_df = CSV.read(reptile_mass_csv, DataFrame)
s_reptile_mass = leftjoin(s, reptile_mass_df; on=:GBIFSpecies, matchmissing=:notequal, makeunique=true, order=:left)
filter(r -> contains(r.Order, "Squa"), reptile_mass_df)

# frugivores_path = "/home/raf/PhD/Mascarenes/Tables/Dryad frugivore occurrence database 1-3-17.xlsx"
# frugivores_xl = XLSX.readxlsx(frugivores_path)
# frugivores_df = DataFrame(XLSX.eachtablerow(frugivores_xl["Frugivore occurrence and traits"]))
frugivores_path = "Dryad frugivore occurrence database 1-3-17.csv"
frugivores_df = CSV.read(frugivores_path, DataFrame)
filter(r -> contains(r.Species_name, "Cylindraspis"), frugivores_df)
s_frugivores = leftjoin(s, frugivores_df; on=:GBIFSpecies, matchmissing=:notequal, makeunique=true, order=:left)
s_frugivores_mass = DataFrames.combine(groupby(s_frugivores, :GBIFSpecies), :Body_mass => mean)
# s_not_frugivores = antijoin(s, frugivores_df; on=:GBIFSpecies=>:Species_name, matchmissing=:notequal, makeunique=true)
# subset(s_not_frugivores, :Origin => ByRow(==("Endemic")))
# subset(s, :Origin => ByRow(==("Endemic")))

# frugivores_df.GBIFSpecies = copy(frugivores_df.Species_name)
# for i in eachindex(frugivores_df.Species_name) 
#     sp = frugivores_df.Species_name[i]
#     ismissing(sp) && continue
#     @show sp
#     match = species_match(sp)
#     isnothing(match) && continue
#     frugivores_df.GBIFSpecies[i] = match.species 
# end
# frugivores_df.GBIFSpecies
# CSV.write(frugivores_path, frugivores_df)

# subset(s_not_frugivores, :Origin => ByRow(==("Endemic")))

sort!(s_elton, :GBIFSpecies)
sort!(s_avonet, :GBIFSpecies)
sort!(s_pantheria, :GBIFSpecies)
sort!(s_frugivores_mass, :GBIFSpecies)
sort!(s_bird_mass, :GBIFSpecies)
sort!(s_reptile_mass, :GBIFSpecies)
sort!(s, :GBIFSpecies)
# filter(r -> ismissing(r.Mass), s)

function combine_mass(a, b)
    @assert length(a) == length(b)
    map(a, b) do s1, s2
        if ismissing(s1) 
            ismissing(s2) ? missing : s2
        else
            s1
        end
    end
end
s.Mass .= missing
s.Mass = combine_mass(s.Mass, s_avonet.Mass_1)
s.Mass = combine_mass(s.Mass, s_elton.BodyMass_Value)
s.Mass = combine_mass(s.Mass, s_pantheria.AdultBodyMass_g)
s.Mass = combine_mass(s.Mass, s_reptile_mass[!, "Weight (g)"])
s.Mass = combine_mass(s.Mass, s_frugivores_mass.Body_mass_mean)
s.Mass = combine_mass(s.Mass, s_bird_mass.filledmass)
subset(s, :Mass => x -> ismissing.(x), :Origin => ByRow(==("Endemic")))

s.Mass_Avonet = s_avonet.Mass_1 
s.Mass_Pantheria = s_pantheria.AdultBodyMass_g 
s.Mass_Elton = s_elton.BodyMass_Value 
s.Mass_Frugivores = s_frugivores_mass.Body_mass_mean
s.Mass_bird_mass = s_bird_mass.filledmass
s.Mass_reptile_mass = s_reptile_mass[!, "Weight (g)"]
s.Mass

masses = (s.Mass_Avonet, s.Mass_Pantheria, s.Mass_Elton, s.Mass_Frugivores, s.Mass_bird_mass, s.Mass_reptile_mass)
collect(skipmissing(s.Mass_Frugivores))
map(masses) do m
    8196 in m
end

# filter(r -> ismissing(r.Mass) && r.Origin == "Endemic", s) |> pager

# elton_end = filter(r -> r.Origin == "Endemic" && !ismissing(r.BodyMass_Value), s_elton)
# elton_end = filter(r -> r.Origin == "Endemic" && !ismissing(r.Mass), s)
# elton_end = filter(r -> r.Origin == "Endemic", s)
# elton_end = filter(r -> r.Origin == "Endemic" && ismissing(r.Mass), s)
# filter(avonet.Species1) do sp
#     occursin("Nycticorax", sp)
# end


combine_csv = "/home/raf/Data/Traits/Combine/COMBINE_archives/trait_data_imputed.csv"
combine = CSV.read(combine_csv, DataFrame; normalizenames=true)

combine_csv = "/home/raf/Data/Traits/Combine/COMBINE_archives/taxonomy_crosswalk.csv"
combine = CSV.read(combine_csv, DataFrame; normalizenames=true)

catoflife_tsv = "/home/raf/Data/Taxonomy/CatalogueOfLife/Taxon.tsv"
catoflife = CSV.read(catoflife_tsv, DataFrame; normalizenames=true, missingstring=nothing)
names(catoflife)
catoflife.Binomial

full_names = map(catoflife.dwc_genericName, catoflife.dwc_specificEpithet, catoflife.dwc_infraspecificEpithet) do g, s, i
    strip("$g $s $i")
end
catoflife.Scientific = full_names
length(catoflife.Scientific)
length(union(catoflife.Scientific))

syn = leftjoin(s, catoflife; on=:Species=>:Scientific)
names(catoflife)
for i in 1:nrow(s)
    fam = species_search(s.Species[i]).family
    length(fam) > 0 || continue
    if ismissing(s.Family[i])
        f = first(skipmissing(fam))
        @show s.Family[i] f
        s.Family[i] = f
    # elseif !ismissing(sp[1].family)
        # s.Family[i] == sp[1].family || @show s.Family[i] sp[1].family
    end
end
s.Family

length(collect(skipmissing(syn.dwc_specificEpithet)))
length(collect(skipmissing(syn.dwc_specificEpithet)))

# leftjoin(s, pantheria; on=:Species=>:MSW05_Binomial) |> pager
# leftjoin(s, avonet; on=:Species=>:Species1, makeunique=true) |> pager
# leftjoin(s, lizzard; on=:Species=>:Binomial, makeunique=true, matchmissing=:notequal) |> pager

s.has_trait_data .|= in.(s.Species, Ref(pantheria.MSW05_Binomial))
s.has_trait_data .|= in.(s.Species, Ref(avonet.Species1))
s.has_trait_data .|= in.(s.Species, Ref(lizzard.Binomial))
s.has_trait_data .|= in.(s.Species, Ref(elton_bird.Scientific)) # only adds 4
s.has_trait_data .|= in.(s.Species, Ref(combine.iucn2020_binomial)) # only adds 1
# s.has_trait_data .|= in.(s.Species, Ref(elton_mammal.Scientific)) # Adds nothing 
sum(skipmissing(s.has_trait_data))
# filter(r -> ismissing(r.has_trait_data), s) |> pager
#

# filter(e) do r
    # ismissing(r.rod) && !ismissing(r.rod_extinct)
# end

# species_params = map(s.Rmax, s.Max_Density, s.mus_extinct, s.reu_extinct, s.rod_extinct) do rmax, max_density, mus_extinct, reu_extinct, rod_extinct
#     (; rmax, max_density, mus_extinct, reu_extinct, rod_extinct)
# end |> Tuple |> NamedVector{get_species_names(s)}

# island_params = map(island_extinct_names, island_names) do keys, island
#     spec = species_params[keys]
#     map(spec) do s
#         extinct = s[Symbol(island, "_extinct")]
#         base = (; s[(:rmax, :max_density, :hunting_suscept, :cat_suscept, :rat_suscept)]..., extinct)
#     end
# end

# island_populations = map(island_params) do params
#     v = NamedVector(map(_ -> 100.0f0, params))
#     fill(v, 100, 100)
# end

# island_columns = map(island_params) do params
#     k = keys(first(params))
#     map(k) do key
#         NamedVector(map(r -> r[key], params))
#     end |> NamedTuple{k}
# end

# island_columns.mus.max_density
# island_params.mus


# Population Based??
# growth_rules = map(island_columns) do island
#     growth = Dispersal.LogisticGrowth{:populations,:populations}(;
#         rate=island.rmax,
#         carrycap=island.max_density,
#         timestep=1,
#         nsteps_type=Float32,
#     )
# end

# hunting = let hunting_suscept = hunting_suscept
#     Cell{:populations,:populations}() do data, pops, I
#         hp = get(data, Aux{:hunting_pressure}(), I)
#         pops .* (1 .- hunting_suscept .* hunting_pressure)
#     end
# end

# cat_predation = let cat_suscept = cat_suscept, 
#                     weights = weights
#     Cell{Tuple{:populations,:cats}}() do data, (pops, cats), I
#         max_predation = pops .* cat_suscept .* cats
#         max_nutrition = weights .* meat_kj .* potential_predation
#         max_required = cats * cat_K .* max_nutrition
#         predation = max_predation .* frac
#         cats = 
#         newpops .- predation, 
#     end
# end

# ruleset = Ruleset(growth, spread, cat_predation, rat_predation, hunting)

# Presence/absense Based??

# settlement = (
#     1598, "Port de Warwick/Grand port used as stopover", "Dutch", -20.3754, 57.7228,
#     1606, "Port Louie used as stopover", "Dutch", -20.1597, 57.5012,
#     1635, "First garrison at Port de Warwick", "Dutch", -20.3754, 57.7228, "Ebony clearing in this time",
#     1645, "Road build in Flacq",
#     1655, "Three settlements: Grand port bay, Flacq, Trou dEau Douce",
#     ~1700, "Flacq and Black River settlements",
#     1721, "Colony at port louie founded", "French", -20.1597, 57.5012,
#     1806, "Mahebourge founded", "French", -20.4045, 57.7028,
#     1814, "Mauritius ceded to britain", "English",
# )

