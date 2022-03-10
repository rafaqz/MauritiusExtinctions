using CSV
using DataFrames
using Plots

# Population
human_pop = CSV.File(joinpath(workdir, "Data/Population/Population.csv")) |> DataFrame
sugar_cane = CSV.File(joinpath(workdir, "Data/Population/Sugarcane.csv")) |> DataFrame
human_pop.Population .*= 1000
# Plots.plot(human_pop.Year, human_pop.Population)
# Plots.plot(sugar_cane.Year, sugar_cane.Area)


# Species distributions
#
# GBIF alternatives
# Arctos, neotoma, vertnet
using GBIF, Plots, Rasters, IntervalSets
species = CSV.File("/home/raf/PhD/Mauritius/mascarine_species.csv") |> DataFrame

# Google search all species, opening the next when the browser is closed.
for row in eachrow(subset(species, :Origin => ByRow(==("Endemic"))))
    ismissing(row["Common name"]) || ismissing(row["Species"]) && continue
    placenamed = false
    for place in ("RéuReunion", "Mauritius", "Rodrigues", "Mascarene")
        placenamed = occursin(place, row["Common name"]) 
        @show row["Common name"] place placenamed
        placenamed && break
    end
    placenamed && continue
    search = replace(row["Species"], " " => "+")
    run(`chromium https\://www.google.com/search\?q=$(search)`)
    # run(`chromium https\://www.google.com/search\?tbm=isch\&q=$(search)`)
end

for sp in species.Species
    ismissing(sp) && continue
    search = replace(sp, " " => "+")
    run(`chromium https\://www.google.com/search\?q=$(search)`)
    # run(`chromium https\://www.google.com/search\?tbm=isch\&q=$(search)`)
end

# Check in Hume
humepath = "/home/raf/zotero_library/Journal Article/undefined/2013/Hume_2013_A*"
for sp in species.Species
    ismissing(sp) && continue
    run(`evince --find=$(sp) "/home/raf/zotero_library/Journal Article/undefined/2013/Hume_2013_A synopsis of the pre-human avifauna of the Mascarene Islands.pdf"`)
end


endemics = DataFrames.subset(species, :Origin => ByRow(==("Endemic")); skipmissing=true)
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

while length(obs) < size(obs)
    @show length(obs) size(obs)
    occurrences!(obs)
end
