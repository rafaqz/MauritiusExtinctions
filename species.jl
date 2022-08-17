using XLSX, DataFrames
using Chain
using Plots, StatsPlots
using DimensionalData, DimensionalData.LookupArrays
using OrderedCollections
using IntervalSets
using DimensionalData
using DimensionalData.LookupArrays

pyplot()
include("tabular_data.jl")
include("plots.jl")
include("lostland_filter.jl")

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
mascarine_species = CSV.File(joinpath(workdir, "mascarine_species.csv")) |> DataFrame
xlfile = joinpath(datadir, "LostLand/Mauritius_Lost Land of the Dodo_tables_translated symbols.xlsx")
# run(`libreoffice $xlfile`)
xl = XLSX.readxlsx(xlfile)
pops = map(lost_land_appendices) do island
    map(island) do sheetname
        @show sheetname
        filter_population(as_dataframe(xl, sheetname))
    end
end;

timelineplots = map(NamedTuple{keys(lost_land_appendices)}(keys(lost_land_appendices))) do island
    map(NamedTuple{keys(first(lost_land_appendices))}(keys(first(lost_land_appendices)))) do category
        p = plottimeline(mascarine_species, island, category)
        savefig(p, "$(island)_$category.png")
        p
    end
end
introplot = plot(plant_introductions;
    legend=:topleft,
    link=:x,
)
plot(timelineplots.reu..., introplot;
    layout=grid(3, 1, heights=(0.45,0.45,0.1)),
    size=(2000, 2000),
    link=:x,
)

species = map(pops) do island
    dfn = island.native
    extinction = map(names(dfn)) do name
        obs = dfn[!, name]
        i = findlast(x -> (!ismissing(x) && x != 0), obs)
        name => (isnothing(i) ? missing : dims(obs, Ti())[i])
    end |> OrderedDict
    dfa = island.alien
    introduction = map(names(dfa)) do name
        obs = dfa[!, name]
        i = findfirst(x -> !ismissing(x), obs)
        name => (isnothing(i) ? missing : dims(obs, Ti())[i])
    end |> OrderedDict
    (; extinction, introduction)
end

species_years = map(species) do island
    map(island) do class
        filter(sort(collect(class); by=last)) do x
            !ismissing(x[2]) && x[2] < 2000
        end
    end
end

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

records = map(endemics[!, :Species]) do sp
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
