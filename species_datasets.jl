using StaticArrays
using Statistics
using StatsBase
using Shapefile
using CSV, DataFrames, XLSX, TerminalPager
using GBIF2
using CairoMakie
using GLMakie
using Loess
using Random

function set_gbif_species!(df, specieskey)
    if !("GBIFSpecies" in names(df))
        df.GBIFSpecies .= ""
    end
    specvec = collect(getproperty(df, specieskey))
    for i in eachindex(specvec)
        current = df.GBIFSpecies[i]
        ismissing(current) || current == "" || continue
        sp = specvec[i]
        ismissing(sp) && continue
        match = GBIF2.species_match(sp)
        df.GBIFSpecies[i] = if isnothing(match) || ismissing(match.species)
            specvec[i]
        else
            match.species
        end
    end
end

function meanmass(xs)
    xs1 = filter(skipmissing(xs)) do x
        x > 0
    end
    if length(xs1) > 0
        mean(xs1)
    else
        missing
    end
end

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

# names(s)
# set_gbif_species!(s, :Species)
# CSV.write(hawaii_csv_path, s)

extinctions_csv_path = "/home/raf/PhD/Mascarenes/MauritiusExtinctions/tables/IUCN_extinctions.csv"
extinctions_download = "/home/raf/Downloads/IUCN Extinctions - assessments_gbif.csv"
isfile(extinctions_download) && mv(extinctions_download, extinctions_csv_path; force=true)
s = CSV.read(extinctions_csv_path, DataFrame; types=Dict(:LocationColonised=>Int, :ArchipelagoColonised=>Int)) |>
    x -> filter(x) do row
    !ismissing(row.phylumName) &&
    row.scientificName != "Chenonetta finschi" && # Extinct from Mauris, not europeans
    # row.kingdomName == "ANIMALIA" &&
    row.className in ("REPTILIA", "AMPHIBIA", "MAMMALIA", "AVES") && # No fish or molluscs
    row.systems != "Marine" && # No marine species like seals or whales
    true
end
set_gbif_species!(s, :scientificName)

# mascarene_species_csv = "/home/raf/PhD/Mascarenes/MauritiusExtinctions/tables/mascarene_species.csv"
# s = CSV.read(mascarene_species_csv, DataFrame)# |>
# @async run(`libreoffice $mascarene_species_csv`)
    # x -> subset(x, :Species => ByRow(!ismissing), :GBIFSpecies => ByRow(!ismissing))
# CSV.write(mascarene_species_csv, s)
#
# iucn_extinct = CSV.read("/home/raf/Data/Extinction/redlist/extinct_species.csv", DataFrame)
# iucn_extinct.threats |> pager
# iucn_csv = "/home/raf/Data/Extinction/redlist/redlist_species_data_39da78ce-d594-4968-8043-489f2765d687/assessments_gbif.csv"
# iucn_df = CSV.read(iucn_csv, DataFrame)
# set_gbif_species!(iucn_df, :scientificName)
# CSV.write(iucn_csv, iucn_df)
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

# lizzard = CSV.read(lizzard_csv, DataFrame; normalizenames=true)
# mass_cols = ["Binomial", "mass_equation_Feldman_et_al_2016_unless_stated_", "intercept", "slope"]
# lizzard_mass = lizzard[!, mass_cols]
# s_lizzard = leftjoin(s, lizzard_mass; on=:GBIFSpecies => :Binomial, matchmissing=:notequal, makeunique=true)

# Define all trait dataframes, with key column names
trait_csvs = (;
    hawaii=(csv="/home/raf/PhD/Mascarenes/MauritiusExtinctions/tables/FE.Case.Tarwater_2020.csv", mass=Symbol("Body.mass_grams"), binomial=:Species),
    mascarene=(csv="tables/mascarene_species.csv", mass=:Mass, binomial=:Species),
    pantheria=(csv="/home/raf/Data/Traits/PanTHERIA/ECOL_90_184/PanTHERIA_1-0_WR05_Aug2008_gbif.csv", mass=:AdultBodyMass_g, binomial=:MSW05_Binomial),
    avonet = (csv="/home/raf/Data/Traits/Avonet/ELEData/ELEData/TraitData/AVONET1_BirdLife_gbif.csv", mass=:Mass, binomial=:Species1),
    # lizzard = (csv="/home/raf/Data/Traits/Lizards/Appendix S1 - Lizard data version 1.0.csv" binomial=:XX),
    elton_mammal = (csv="/home/raf/Data/Traits/EltonTraits/MamFuncDat_gbif.csv", mass=:BodyMass_Value, binomial=:Scientific),
    elton_bird = (csv="/home/raf/Data/Traits/EltonTraits/BirdFuncDat_gbif.txt", mass=:BodyMass_Value, binomial=:Scientific),
    reptile_mass = (csv="/home/raf/PhD/Mascarenes/Tables/Reptile body mass database Meiri 2010_gbif.csv", mass=Symbol("Weight (g)"), binomial=:Name),
    bird_mass = (csv="/home/raf/PhD/Mascarenes/Tables/Bird Mass filled (Jan 22 2015)_WDK_gbif.csv", mass=:filledmass, binomial=:BirdLife_SpecName),
    frugivores = (csv="tables/Dryad frugivore occurrence database 1-3-17.csv", mass=:Body_mass, binomial=:Species_name),
)

heinen_csv = "/home/raf/PhD/Mascarenes/MauritiusExtinctions/tables/Heinen_extinct_terrestrial_vertebrates.csv"
heinen_mass = CSV.read(heinen_csv, DataFrame;
    missingstring="NA", types=Dict(:Mean_Body_Mass_Heinen_gram => Float64)
) |> x -> select(x, [:GBIFSpecies, :Mean_Body_Mass_Heinen_gram])

# Open all csvs as DataFrames
trait_dfs = map(trait_csvs) do props
    df = CSV.read(props.csv, DataFrame)
    merge((; df), props)
end

# Make sure to update all the names from GBIF
map(trait_dfs, keys(trait_dfs)) do props, key
    set_gbif_species!(props.df, props.binomial)
end

# Combined all the trait dataframes into one
mega_mass_df = map(trait_dfs, keys(trait_csvs)) do (; df, mass, binomial), db
    df.Database .= db
    # Standardize mass column
    dropmissing(select(df, :GBIFSpecies, binomial => :OriginalBinomial, mass => :Mass, :Database))
end |> splat(vcat) |> sort

# Split and combine taking the mean, will tracking the source databases
mean_mass_df = DataFrames.combine(
    groupby(mega_mass_df, :GBIFSpecies),
    [:Mass => meanmass => :Mass_mean, :Database => Tuple => :Mass_sources],
) |> sort


# Split out the Genus field to another column
mean_mass_df.Genus = map(mean_mass_df.GBIFSpecies) do s
    ismissing(s) ? missing : split(s, ' ')[1]
end
s.Genus = map(s.GBIFSpecies) do s
    ismissing(s) ? missing : split(s, ' ')[1]
end

genus_mean_mass_df = DataFrames.combine(
    groupby(mean_mass_df, :Genus),
    [:Mass_mean => meanmass => :Genus_mass_mean, :Mass_sources => Tuple => :Genus_mass_sources],
) |> sort

s_mass = leftjoin(s, mean_mass_df; on=:GBIFSpecies, matchmissing=:notequal, makeunique=true) |>
     x -> leftjoin(x, genus_mean_mass_df; on=:Genus, matchmissing=:notequal) |>
     x -> leftjoin(x, heinen_mass; on=:GBIFSpecies)

s_mass.EstimatedMass = map(s_mass.Mass_mean, s_mass.Genus_mass_mean, s_mass.Mean_Body_Mass_Heinen_gram) do mm, gm, hm
    # First check dataset species mean
    x = if ismissing(mm) || isnan(mm)
        # Then Heinen, otherwise use genus mean
        ismissing(hm) ? gm : hm
    else
        mm
    end
    (ismissing(x) || isnan(x)) ? missing : x
end;
s_mass.EstimatedMass |> skipmissing |> collect |> length
s_mass.colonised = map(s_mass.ArchipelagoColonised, s_mass.LocationColonised) do a, i
    ismissing(i) ? a : i
end
s_mass.isisland = s_mass.Island .== "Yes"
s_mass.wasuninhabited = map(s_mass.ArchipelagoPreviouslyInhabited .== ("No",), s_mass.LocationPreviouslyInhabited .== ("No",)) do a, i
    if ismissing(a)
        if ismissing(i)
            false
        else
            i
        end
    elseif ismissing(i)
        a
    else
        a || i
    end
end

s_mass
# simplified = select(s_mass, [:scientificName, :GBIFSpecies, :Archipelago, :Location, :Mass_mean, :Genus_mass_mean, :Island])

weigelt_csv = "/home/raf/Data/Extinction/Islands/Weigelt/Weigelt_etal_2013_PNAS_islanddata.csv"
# run(`libreoffice $weigelt_csv`)
weigelt_islands = CSV.read(weigelt_csv, DataFrame)
names(weigelt_islands)

s_colonised = dropmissing(s_mass, :colonised)
# s_weigelt = leftjoin(s_colonised, weigelt_islands; on=:WeigeltID=>:ID, matchmissing=:notequal, makeunique=true)
s1 = dropmissing(s_mass, [:EstimatedMass, :yearLastSeen_cleaned, :Area])
s_no_mass = subset(s_colonised, :EstimatedMass => x -> ismissing.(x))
sort(collect(s_no_mass.Location |> countmap); by=last)

sort(s1.Location |> countmap |> pairs |> collect; by=last)
sort(s1.Archipelago |> countmap |> pairs |> collect; by=last)
sort(s1.SuperArchipelago |> countmap |> pairs |> collect; by=last)

function plot_extinctions!(ax, df;
    colonised=nothing,
    names=nothing,
    colordata=:colonised,
    colormap=:viridis
)
    longlived = subset(df, :EstimatedMass => ByRow(>(1e4)))
    xlims!(ax, (1480, 2020))
    ylims!(ax, (1.0, 1e6))
    xs, ys = df.yearLastSeen_cleaned, df.EstimatedMass
    # xs, ys = s1.yearLastSeen_cleaned .- s1.colonised, log.(s1.Mass)
    # ys = shuffle(ys)
    p = Makie.plot!(ax, xs, ys;
        label="Extinctions",
        color=getproperty(df, colordata),
        colorrange=(1500, 1900),
        colormap,
        inspector_label=(_, i, _) -> "$(df.GBIFSpecies[i])\nIsland: $(df.Location[i])\nArea: $(df.Area[i])\nMass: $(ys[i])\nExtinct: $(xs[i])",
    )
    if names == :text
        Makie.text!(ax, xs, ys; text=df.GBIFSpecies)
    end
    # p = Makie.plot!(ax, longlived.yearLastSeen_cleaned, longlived.EstimatedMass; color=(:yellow, 0.5), label="Long lived")#; color=df.colonised, colormap=:viridis)

    model = loess(xs, log.(ys), span=1.0, degree=2)
    us = range(extrema(xs)...; step=5)
    vs = exp.(predict(model, us))
    Makie.lines!(ax, us, vs; label="Loess regression", color=(:black, 0.8))
    Makie.hlines!(ax, median(ys); label="Median mass")
    if !isnothing(colonised)
        Makie.vlines!(ax, df.colonised; color=colonised, label="Colonisation")
    end
    DataInspector(ax)
    return nothing
end

sort(collect(s1.Archipelago |> countmap); by=last)

not_mauris = :GBIFSpecies => ByRow(!=("Chenonetta finschi"))
subset_queries = (;
    all=(title="All colonised", query=()),
    islands=(title="All Islands", query=(:isisland,)),
    continents=(title="Continents", query=(:isisland => .!,)),
    inhabited_islands=(title="Inhabited", query=(not_mauris, :isisland, :wasuninhabited => .!,)),
    uninhabited_islands=(title="Uninhabited", query=(:isisland, :wasuninhabited,)),
    indian_ocean=(title="Indian ocean", query=(:SuperArchipelago=>ByRow(==("Indian Ocean")),)),
    mascarenes=(title="Mascarenes", query=(:Archipelago=>ByRow(==("Mascarenes")),)),
    non_mascarene_uninhabited=(title="Non-Mascarene Uninhabited", query=(:isisland, :wasuninhabited, :Archipelago=>ByRow(!=("Mascarenes")),)),
    mauritius=(title="Mauritius", query=(:Location=>ByRow(==("Mauritius")),)),
    reunion=(title="Reunion", query=(:Location=>ByRow(==("Reunion")),)),
    rodrigues=(title="Rodrigues", query=(:Location=>ByRow(==("Rodrigues")),)),
    australia=(title="Australia", query=(:SuperArchipelago=>ByRow(==("Australia")),)),
    new_zealand=(title="New Zealand", query=(:SuperArchipelago=>ByRow(==("New Zealand")),)),
    st_helena=(title="St Helena", query=(:SuperArchipelago=>ByRow(==("St Helena")),)),
    west_indies=(title="West Indies", query=(:SuperArchipelago=>ByRow(==("West Indies")),)),
    hawaiian_islands=(title="Hawaiian Islands", query=(:Archipelago=>ByRow(==("Hawaiian Islands")),)),
    polynesia=(title="Polynesia", query=(:SuperArchipelago=>ByRow(==("Polynesia")),)),
    micronesia=(title="Micronesia", query=(:SuperArchipelago=>ByRow(==("Micronesia")),)),
    galapagos=(title="Galapagos", query=(:Archipelago=>ByRow(==("Galapagos")),)),
    australian_continent=(title="Australian Continent", query=(:isisland => .!, :SuperArchipelago=>ByRow(==("Australia")),)),
    australian_islands=(title="Australian Continent", query=(:isisland, :SuperArchipelago=>ByRow(==("Australia")),)),
    australian_inhabited_islands=(title="Australian Inhabited Islands", query=(:isisland, :wasuninhabited => .!, :SuperArchipelago=>ByRow(==("Australia")),)),
    australian_uninhabited_islands=(title="Australian Uninhabited Islands", query=(:isisland, :wasuninhabited, :SuperArchipelago=>ByRow(==("Australia")),)),
    all_early=(title="All Early Colonisation", query=(:isisland, :colonised=>ByRow(<(1750)),)),
    all_late=(title="All Late Colonisation", query=(:isisland, :colonised=>ByRow(>=(1750)),)),
    inhabited_early=(title="Inhabited Early Colonisation", query=(not_mauris, :isisland, :wasuninhabited => .!, :colonised=>ByRow(<(1750)),)),
    inhabited_late=(title="Inhabited Late Colonisation", query=(:isisland, :wasuninhabited => .!, :colonised=>ByRow(>=(1750)),)),
    uninhabited_early=(title="Uninhabited Early Colonisation", query=(:isisland, :wasuninhabited, :colonised=>ByRow(<(1750)),)),
    uninhabited_late=(title="Uninhabited Late Colonisation", query=(:isisland, :wasuninhabited, :colonised=>ByRow(>=(1750)),)),
)
subs = map(subset_queries) do qs
    df = subset(s1, qs.query...; skipmissing=true)
    merge(qs, (; df))
end

subs.uninhabited_islands.df.Location |> union
subs.uninhabited_islands.df.SuperArchipelago |> union
subs.uninhabited_islands.df.SuperArchipelago |> countmap
subs.uninhabited_islands.df.Archipelago |> countmap
subs.uninhabited_islands.df |> nrow

subset_layout = [
    :islands                 :all_early         :all_late #nothing
    :inhabited_islands   :inhabited_early   :inhabited_late #:west_indies
    :uninhabited_islands :uninhabited_early :uninhabited_late #nothing
]
kw = (; yscale=log10, xlabel="Year last seen", ylabel="Mass")
fig = Figure(; size=(1600,900));
colormap = :managua
colordata = :Area
axs = map(subset_layout, CartesianIndices(subset_layout)) do key, I
    if isnothing(key)
        ax = Axis(fig[Tuple(I)...])
        xlims!(ax, (1400, 2020))
        ylims!(ax, (1.0, 1e6))
    else
        @show key
        subset = subs[key]
        ax = Axis(fig[Tuple(I)...]; title=subset.title, kw...)
        I[1] == size(subset_layout, 1) || hidexdecorations!(ax; grid=false)
        I[2] == 1 || hideydecorations!(ax; grid=false)
        plot_extinctions!(ax, subset.df;
            # colonised=(:black, 0.02),
            names=:tooltip,
            colordata, colormap
        )
    end
    ax
end
linkaxes!(axs...)
Makie.Colorbar(fig[UnitRange(axes(subset_layout, 1)), size(subset_layout, 2) + 1];
    colormap,
    colorrange=extrema(getproperty(s1, colordata)),
    label="Colonisation date",
)
# axislegend(axs[1]; position=:lt)
fig[0, :] = Label(fig, "Patterns of mass, extinction date, and human habitation", fontsize=20)
display(fig)
# save("images/dome_extinction.png", fig)

# using GLMakie
GLMakie.activate!()
fig = Figure(; size=(900,600));
kw = (; yscale=log10, xlabel="Year of extinction", ylabel="Mass")
ax1 = Axis(fig[1, 1]; kw...)
plot_extinctions!(ax1, subs.australian_continent.df; names=true)
display(fig)

# plot_extinctions!(ax4, subs.rodrigues; colonised=true)


# Makie.vlines!([1890, 190, 1910]; color=:gray)
# Makie.vlines!([1950]; color=:red)
# Makie.vlines!([1955]; color=:blue)


# g = groupby(s1, :Location; skipmissing=true)
# keys(g)
# sort(keys(g) .=> [nrow(x) for x in g]; by=last)
# g[("Norfolk Island Group",)]
# s1.Archipelago |> union |> sort

# s1 = subset(dropmissing(s, [:Mass, :yearLastSeen_cleaned, :colonised]), :Island => ByRow(==("Yes")))


# Late introductions
# Mongoose introduced 1870, 1900, 1910
# Brown tree snake introduced to guam 1945-1952
# Rosy wolf snail 1955

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

