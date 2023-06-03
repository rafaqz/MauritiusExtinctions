using DynamicGrids, Dispersal, LandscapeChange, ModelParameters
using StaticArrays 
using Shapefile
using CSV, DataFrames, XLSX, TerminalPager
using GBIF2

includet("raster_common.jl")

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

mascarine_species_csv = "/home/raf/PhD/Mascarenes/Tables/mascarine_species.csv"
s = CSV.read(mascarine_species_csv, DataFrame) |> 
    x -> subset(x, :Species => ByRow(!ismissing))

# s.GBIFSpecies = copy(s.Species)
# for i in eachindex(s.Species) 
#     sp = s.Species[i]
#     ismissing(sp) && continue
#     match = species_match(sp)
#     isnothing(match) && continue
#     s.GBIFSpecies[i] = match.species 
# end
# s.GBIFSpecies

# CSV.write(mascarine_species_csv, s)
iucn_extinct = CSV.read("/home/raf/Data/Extinction/redlist/extinct_species.csv", DataFrame)
iucn2 = CSV.read("/home/raf/Data/Extinction/redlist/redlist_species_data_39da78ce-d594-4968-8043-489f2765d687/assessments.csv", DataFrame)

names(iucn2)

joined = leftjoin(s, iucn2; on=:GBIFSpecies=>:scientificName, matchmissing=:notequal)
joined[!, [:Species, :threats]] |> pager
broadcast(string, joined.Species, Ref(": "), joined.threats) |> pager


sp_nothing = filter(x -> isnothing(x[2]), sp_pairs)
sp = filter(x -> !isnothing(x[2]), sp_pairs)
gbif_sp = DataFrame(last.(sp)
first.(sp) .=> gbif_sp.species
first.(sp)[first.(sp) .!== gbif_sp.species]

redlist_extinct_csv = "/home/raf/Data/Traits/redlist_species_data_1f74a1f8-0b29-4567-9766-046807e966ca/taxonomy.csv"
redlist = CSV.read(redlist_extinct_csv, DataFrame; normalizenames=true)

pantheria_csv = "/home/raf/Data/Traits/PanTHERIA/ECOL_90_184/PanTHERIA_1-0_WR05_Aug2008.txt"
pantheria = CSV.read(pantheria_csv, DataFrame; normalizenames=true, quoted=false)
pantheria_mass = pantheria[!, [:MSW05_Binomial, :AdultBodyMass_g]]
s_pantheria = leftjoin(s, pantheria_mass; on=:GBIFSpecies=>:MSW05_Binomial, matchmissing=:notequal, makeunique=true)
names(pantheria)
sort!(s_pantheria, :GBIFSpecies)

avonet_csv = "/home/raf/Data/Traits/Avonet/ELEData/ELEData/TraitData/AVONET1_BirdLife.csv"
avonet = CSV.read(avonet_csv, DataFrame; normalizenames=true)
avonet_mass = avonet[!, [:Species1, :Mass]]
s_avonet = leftjoin(s, avonet_mass; on=:GBIFSpecies=>:Species1, matchmissing=:notequal, makeunique=true)
names(avonet)

lizzard_csv = "/home/raf/Data/Traits/Lizards/Appendix S1 - Lizard data version 1.0.csv"
lizzard = CSV.read(lizzard_csv, DataFrame; normalizenames=true)
lizzard_mass = avonet[!, [:Species1, :Mass]]
names(lizzard)

mass_cols = ["Binomial", "mass_equation_Feldman_et_al_2016_unless_stated_", "intercept", "slope"]
lizzard_mass = lizzard[!, mass_cols]
s_lizzard = leftjoin(s, lizzard_mass; on=:GBIFSpecies=>:Binomial, matchmissing=:notequal, makeunique=true)
lizzard_end = filter(r -> !ismissing(r.intercept) && r.Origin == "Endemic", s_lizzard)
# lizzard_end |> pager

elton_bird_csv = "/home/raf/Data/Traits/EltonTraits/BirdFuncDat.txt"
elton_bird = CSV.read(elton_bird_csv, DataFrame; normalizenames=true)
names(avonet)

elton_mammal_csv = "/home/raf/Data/Traits/EltonTraits/MamFuncDat.txt"
elton_mammal = CSV.read(elton_mammal_csv, DataFrame; normalizenames=true)

elton_mass = vcat(elton_mammal[!, [:Scientific, :BodyMass_Value]], elton_bird[!, [:Scientific, :BodyMass_Value]])
s_elton = leftjoin(s, elton_mass; on=:GBIFSpecies=>:Scientific, matchmissing=:notequal)
sort!(s_elton, :GBIFSpecies)
sort!(s_avonet, :GBIFSpecies)
sort!(s_pantheria, :GBIFSpecies)
sort!(s, :GBIFSpecies)
filter(r -> ismissing(r.Mass), s)

function combine_mass(a, b)
    map(a, b) do s1, s2
        if ismissing(s1) 
            ismissing(s2) ? missing : s2
        else
            s1
        end
    end
end

s.Mass = combine_mass(s.Mass, s_avonet.Mass_1)
s.Mass = combine_mass(s.Mass, s_elton.BodyMass_Value)
s.Mass = combine_mass(s.Mass, s_pantheria.AdultBodyMass_g)
s.Mass_Avonet = s_avonet.Mass_1 
s.Mass_Pantheria = s_pantheria.AdultBodyMass_g 
s.Mass_Elton = s_elton.BodyMass_Value 

filter(r -> ismissing(r.Mass) && r.Origin == "Endemic", s) |> pager

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
filter(r -> ismissing(r.has_trait_data), s) |> pager




island_tables = map((mus=:mus, reu=:reu, rod=:rod)) do key
    subset(s, key => x -> .!ismissing.(x))
end
map(length ∘ collect ∘ skipmissing, (s.mus_extinct, s.reu_extinct, s.rod_extinct))
# filter(e) do r
    # ismissing(r.rod) && !ismissing(r.rod_extinct)
# end

get_species_names(table) = Tuple(Symbol.(replace.(skipmissing(table.Species), Ref(" " => "_"))))
island_extinct = map(island_tables, island_names) do table, name
    subset(table, Symbol(name, "_extinct") => ByRow(!ismissing))
end
island_names = NamedTuple{keys(island_species)}(keys(island_species))
island_extinct_names = map(get_species_names, island_extinct)
species_params = map(s.Rmax, s.Max_Density, s.mus_extinct, s.reu_extinct, s.rod_extinct) do rmax, max_density, mus_extinct, reu_extinct, rod_extinct
    (; rmax, max_density, mus_extinct, reu_extinct, rod_extinct, hunting_suscept=1.0, cat_suscept=1.0, rat_suscept=1.0)
end |> Tuple |> NamedVector{get_species_names(s)}

island_params = map(island_extinct_names, island_names) do keys, island
    spec = species_params[keys]
    map(spec) do s
        extinct = s[Symbol(island, "_extinct")]
        base = (; s[(:rmax, :max_density, :hunting_suscept, :cat_suscept, :rat_suscept)]..., extinct)
    end
end

island_populations = map(island_params) do params
    v = NamedVector(map(_ -> 100.0f0, params))
    fill(v, 100, 100)
end

island_columns = map(island_params) do params
    k = keys(first(params))
    map(k) do key
        NamedVector(map(r -> r[key], params))
    end |> NamedTuple{k}
end

island_columns.mus.max_density
island_params.mus

init_pops = map(island_columns, masks) do params, mask
    fill(map(Float32, params.max_density), size(mask)) .* mask
end


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

using DynamicGrids, Dispersal, LandscapeChange, ModelParameters
using StaticArrays 
using CSV, DataFrames, XLSX, TerminalPager
using Rasters

includet("raster_common.jl")
pred_df = CSV.read("animals.csv", DataFrame)

revmasks = map(masks) do A
    reverse(A; dims=Y)
end
ag_masks = map(revmasks) do mask
    Rasters.aggregate(Rasters.Center(), mask, 10)
end

pred_keys = Tuple(Symbol.(pred_df.name))
pred_names = NamedVector{pred_keys}(pred_keys)
pred_rmax = NamedVector{pred_keys}(Tuple(pred_df.rmax))
pred_max_density = NamedVector{pred_keys}(Tuple(pred_df.rmax))
pred_spread = InwardsDispersal{:pred_pops}(
    radius=1,
    formulation=ExponentialKernel()
)
pred_growth = LogisticGrowth{:pred_pops}(; 
    rate=pred_rmax,
    carrycap=pred_max_density,
    timestep=1
)
hunting_suscept = rand(SVector{10})
habitat_requirement = rand(SVector{10})
pred_pops = map(ag_masks) do m
    map(_ -> map(_->0.0, pred_names), m)
end
pred_suscept = map(hunting_suscept) do hs
    rand()
end

# Every species is everywhere initially, in this dumb model
endemic_presences = map(ag_masks) do mask
    map(mask) do m
        map(_ -> m, hunting_suscept) 
    end
end

risks = let hunting_suscept=hunting_suscept, pred_suscept=pred_suscept
    Cell{Tuple{:pred_pops,:endemic_presences},:endemic_presences}() do data, (pred_pops, endemic_presences), I
        # hp = get(data, Aux{:hunting_pressure}(), I)
        map(endemic_presences, hunting_suscept, pred_suscept) do present, hs, ps
            if present
                # rand() < hp * hs & rand() < sum(map(*, ps, pred_pops))
                rand() < sum(map(*, ps, pred_pops))
                # ... etc 
            else
                false
            end
        end
    end
end

habitat = let habitat_requirement = habitat_requirement
    Cell{:presences}() do data, presences, I
        hp = get(data, Aux{:uncleared}(), I)
        map(presences, habitat_requirement) do present, h
            if present
                rand() < hp * hs
            else
                false
            end
        end
    end
end

tspan = 1600:2020
ruleset = Ruleset(pred_growth, pred_spread, risks)
inits = map(pred_pops, endemic_presences) do pred_pops, endemic_presences
    (; pred_pops, endemic_presences)
end
map(size, inits.mus)
inits.mus.endemic_presences |> parent

output = ResultOutput(inits.mus;
    mask=ag_masks.mus, # aux,
    tspan, 
)
@time sim!(output, ruleset)

output = MakieOutput(init_state;
    aux,
    fps=10,
    tspan,
    store=false,
    mask=ag_masks.mus,
    boundary=Remove(),
    padval=0,
    ruleset,
    sim_kw=(; printframe=true),
) do fig, frame
    axis = Axis(fig[1, 1])
    species = Observable(Array(frame[].landcover))
    i = Observable(1)
    on(frame) do f
        species[] = getindex.(f, i[])
        notify(landcover)
    end
    hm = Makie.image!(axis, landcover; colorrange=(-0.5, 5.5), colormap=:inferno)
    return nothing
end
