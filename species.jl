using CSV
using DataFrames
using XLSX
using OrderedCollections
using Chain
using Plots, StatsPlots
using DimensionalData, DimensionalData.LookupArrays
using IntervalSets
const GI = GeoInterface

exitinct_plant_path = "/home/raf/Data/Extinction/Plants/41559_2019_906_MOESM3_ESM.xlsx"
extinct_plants_xl = XLSX.readxlsx(exitinct_plant_path)
extinct_plant_df = DataFrame(XLSX.eachtablerow(xl["Supplementary_Data_1"]))
localities = extinct_plant_df.Locality
locality_codes = Iterators.flatmap(s -> split(s, ","), skipmissing(localities)) |> union |> sort
# locality_conversion = Dict{String,Any}()
for lc in locality_codes
    if !haskey(locality_conversion, lc)
        println("Enter code for $lc")
        clipboard(lc)
        input = readline()
        if !(input == "")
            locality_conversion[lc] = input
        end
    end
end

# locality_conversion["WAS"] = "USA, Washington"
locality_table = DataFrame(
    :source_name => first.(collect(pairs(locality_conversion))),
    :iso_name => last.(collect(pairs(locality_conversion))),
)

# Not Easter Island and Juan Fernand are in the same file

# locality_table.iso_name = rstrip.(locality_table.iso_name)
locality_table = CSV.read("extinct_plant_locality_table.csv", DataFrame)
sort!(locality_table, [:island, :iso_name])
CSV.write("extinct_plant_locality_table.csv", locality_table)

borders = map(locality_table.iso_name) do region
    map(split(region, "+")) do subregion
        args = String.(split(subregion, ", "))
        try
            border = GeoInterface.convert(GeometryBasics, GADM.get(args...).geom[1])
            # println("Found $name")
            return border
        catch
            println("Could not find \"$subregion\"")
            return missing
        end
    end
end;
polyvecs = map(borders) do geoms
    if ismissing(geoms) || any(ismissing, geoms)
        missing
    else
        map(geoms) do mp
            if GI.trait(mp) isa GI.PolygonTrait 
                [mp]
            else
                collect(GI.getgeom(mp))
            end
        end |> Iterators.flatten |> collect
    end
end |> skipmissing |> collect;
region_multipolygons = GeometryBasics.MultiPolygon.(polyvecs);
region_multipolygons[1]
using GLMakie
poly(region_multipolygons)
# locality_table.geometry = borders;
# missed = filter(row -> any(ismissing, row.geometry), locality_table)
# sort(missed, :iso_name)
# using GeoJSON
# Make GeoJSON write tables...

includet("tabular_data.jl")
includet("plots.jl")
includet("lostland_filter.jl")

years = 1600:2007

# plants_csv = joinpath(workdir, "tassim_reunion_plant_introductions.csv")
# run(`libreoffice $plants_csv`)
# reu_plants = @chain begin
#     CSV.File(plants_csv)
#     DataFrame
#     dropmissing(_, :Introduction_date)
#     sort(_, :Introduction_date)
# end
# reu_plants[!, :Introduction_date] = map(reu_plants.Introduction_date) do d
#     length(d) >= 4 ? parse(Int, d[end-3:end]) : missing
# end
# reu_plants = dropmissing(reu_plants, :Introduction_date)

# counts = zeros(axes(years))
# for (i, y) in enumerate(years)
#     counts[i] = count(==(y), reu_plants.Introduction_date)
# end
# plant_introductions = DimArray(cumsum(counts), Ti(years); name=:num_alien_plants)

lost_land_appendices = (
    mus=(native="Appendix 2", alien="Appendix 5"),
    reu=(native="Appendix 3", alien="Appendix 6"),
    rod=(native="Appendix 4", alien="Appendix 7"),
)

# Turn classifications from the book into population values 0:3,
# representing extinct/not arrived, rare, common, abundant.
#
# Its not clear how to deal with grouped categories like `rats`

using CSV
mascarene_species = CSV.File(joinpath(workdir, "Tables/mascarene_species.csv")) |> DataFrame
xlfile = joinpath(datadir, "LostLand/Mauritius_Lost Land of the Dodo_tables_translated symbols.xlsx")

# run(`libreoffice $xlfile`)
xl = XLSX.readxlsx(xlfile)
pops = map(lost_land_appendices) do island
    map(island) do sheetname
        @show sheetname
        filter_population(as_dataframe(xl, sheetname))
    end
end;

function extinctions(df)
    extinction = map(names(df)) do name
        obs = df[!, name]
        i = findlast(x -> (!ismissing(x) && x != 0), obs)
        period = if isnothing(i) || i == lastindex(obs) 
            missing 
        else
            parse(Int, df.period[i][end-3:end])
        end
        name => period
    end |> OrderedDict
end

species_stats = map(pops) do island
    alien_introduction = map(names(island.alien)) do name
        obs = island.alien[!, name]
        i = findfirst(x -> !ismissing(x), obs)
        name => isnothing(i) ? missing : parse(Int, island.alien.period[i][end-3:end])
    end |> OrderedDict
    (; 
         alien_introduction, 
         alien_extinction = extinctions(island.alien),
         native_extinction = extinctions(island.native),
    )
end

emptyyear(l) = Array{Union{Int,Missing}}(undef, l) .= missing
len = nrow(mascarene_species)

newcols = (
    mus_alien_introduction = emptyyear(len),
    mus_alien_extinction = emptyyear(len),
    mus_native_extinction = emptyyear(len),
    reu_alien_introduction = emptyyear(len),
    reu_alien_extinction = emptyyear(len),
    reu_native_extinction = emptyyear(len),
)

insertcols!(mascarene_species, pairs(newcols)...)

function _fill!(df, colname, dict)
    for (k, v) in pairs(dict)
        df[df.LostLand_name .=== k, colname] = v
    end
end

species_stats.mus.alien_introduction
df.LostLand_name .=== k

_fill!(mascarene_species, :mus_alien_introduction, species_stats.mus.alien_introduction)
_fill!(mascarene_species, :mus_alien_extinction, species_stats.mus.alien_extinction)
_fill!(mascarene_species, :mus_native_extinction, species_stats.mus.native_extinction)
_fill!(mascarene_species, :reu_alien_introduction, species_stats.reu.alien_extinction)
_fill!(mascarene_species, :reu_alien_extinction, species_stats.reu.alien_extinction)
_fill!(mascarene_species, :reu_native_extinction, species_stats.reu.native_extinction)
collect(skipmissing(mascarene_species.reu_native_extinction))
collect(skipmissing(mascarene_species.mus_native_extinction))

timelineplots = map(NamedTuple{keys(lost_land_appendices)}(keys(lost_land_appendices))) do island
    map(NamedTuple{keys(first(lost_land_appendices))}(keys(first(lost_land_appendices)))) do category
        p = plottimeline(mascarine_species, island, category)
        savefig(p, "$(island)_$category.png")
        p
    end
end

pyplot()
# introplot = plot(plant_introductions;
#     legend=:topleft,
#     link=:x,
# )
plot(timelineplots.mus...;#, introplot;
    layout=grid(3, 1, heights=(0.45,0.45,0.1)),
    size=(2200, 2200),
    link=:x,
)

# species = map(pops) do island
#     dfn = island.native
#     extinction = map(names(dfn)) do name
#         obs = dfn[!, name]
#         i = findlast(x -> (!ismissing(x) && x != 0), obs)
#         name => (isnothing(i) ? missing : dims(obs, Ti())[i])
#     end |> OrderedDict
#     dfa = island.alien
#     introduction = map(names(dfa)) do name
#         obs = dfa[!, name]
#         i = findfirst(x -> !ismissing(x), obs)
#         name => (isnothing(i) ? missing : dims(obs, Ti())[i])
#     end |> OrderedDict
#     (; extinction, introduction)
# end

# species_years = map(species) do island
#     map(island) do class
#         filter(sort(collect(class); by=last)) do x
#             !ismissing(x[2]) && x[2] < 2000
#         end
#     end
# end

declen = 10
decades = 1590:decade:1990
groups = map(species_years) do island
    map(island) do class
        groups = []
        for d in decades
            group = [s for s in class if s[2] in d..(d+declen-1)]
            push!(groups, group)
        end
        groups
    end
end

decade_counts = map(groups) do island
    map(class -> length.(class), island)
end
decade_annotations = map(groups) do island
    map(island) do class
        map(class) do x
            names = first.(x)
            join(names, ", ")
        end
    end
end
labels = [string.(keys(groups))...;;]
count_matrxs = map((introduction=:introduction, extinction=:extinction)) do key
    hcat(map(x -> getfield(x, key), values(decade_counts))...)
end

valigns = (mus=:bottom, reu=:center, rod=:top)
anns = map(decade_counts, decade_annotations, valigns) do cs, ts, valign
    offset = 2.5
    shift = valign == :top ? offset : valign == :bottom ? -offset : 0
    map(cs, ts) do c, t
        t = text.(t; pointsize=6, halign=:left)
        tuple.(0, decades .+ shift, t) 
    end
end
ann = anns.mus
kw = (linealpha=0, orientation=:horizontal, yflip=true, labels)
int = groupedbar(decades, count_matrxs.introduction; 
    xlab="Introductions", legend=:none, kw...
)
annotate!(int, anns.mus.introduction)
annotate!(int, anns.reu.introduction)
annotate!(int, anns.rod.introduction)
ext = groupedbar(decades, count_matrxs.extinction;
    xlab="Extinctions", yticks=:none, kw...
)
annotate!(ext, anns.mus.extinction)
annotate!(ext, anns.reu.extinction)
annotate!(ext, anns.rod.extinction)
plot(int, ext; size=(1000, 2000))#; layout=(2, 1))

# Rough ranking of commonness over time
ranked = map(pops) do origins
    map(origins) do origin
        sort(names(origin) .=> sum.(skipmissing.(eachcol(origin))); by=last, rev=true)
    end
end
ranked.reu.alien

key_invasives = [
    "goats",
    "cattle",
    "pigs",
    "cats",
    "rats",
    "Norway Rat",
    "Sa
    hip Rat",
    "Crab-eating Macacque",
]

# Species distributions
#
# GBIF alternatives
# Arctos, neotoma, vertnet
using GBIF2, Plots, Rasters, IntervalSets

endemics = DataFrames.subset(species.mus, :Origin => ByRow(==("Endemic")); skipmissing=true)
lats, lons = (-22.0, -18.0), (55.0, 58.0)
bounds_flags = "decimalLatitude" => lats, "decimalLongitude" => lons

records = map(species.mus[!, :Species]) do sp
    isnothing(sp) || ismissing(sp) && return sp => missing
    taxon = GBIF.taxon(sp)
    return sp => isnothing(taxon) ? missing : DataFrame(GBIF.occurrences(taxon, "limit"=>300, bounds_flags...))
end |> Dict

for (k, v) in records
    ismissing(v) && continue
    df = DataFrame(v)
    if all(ismissing.(df[!, :latitude]))
        records[k] = missing
    else
        @show k length(collect(skipmissing(df.latitude)))
    end
end

df = DataFrame(records["Gallinula chloropus"])
points = collect(zip(skipmissing(df.longitude), skipmissing(df.latitude)))

prec = Raster(WorldClim{Climate}, :prec; month=1, res="30s")
prec = Raster(CHELSA{Climate}, :prec; month=1, res)
mascarines_prec = prec[X = Interval(lons...), Y=Interval(lats...)]
Plots.plot(mascarines_prec)
scatter!(points; opacity=0.5, markershape=:circle)
obs = occurrences(, "limit"=>300)


# Tortoise
# Flat island/Garbriel island had 6000 tortoises. p 205
tortoise_carry_cap = 6000 / fi_area

# Macaque
# How large were historic populations?
# From Sussman and Tattersall
# 25000 - 35000
# Prefere secondary forest!
macaque_1986_pop = 30000
