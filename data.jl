using CSV
using DataFrames
using GeoJSON
using GADM
using Setfield 
using Shapefile
using RasterDataSources
using Rasters
using Rasters.LookupArrays
using Rasters: set, Between, trim, Band
using GeoInterface
using Plots
using Makie, GLMakie
includet("functions.jl")
years = 1638, 1773, 1835, 1872, 1935, "present"

workdir = "/home/raf/PhD/Mauritius/"

soilraster = Raster(joinpath(workdir, "Data/Generated/warpedsoiltypes.tif"))[Band(1)]
Plots.plot(soilraster)
# lakesraster = Raster("warpedlakes.tif")[Band(1)]
# elevationraster = Raster("warpedelevation.tif")[Band(1)]
landuse_snapshots = map(1:6) do i
    Raster("warped_landuse_snapshot_$i.tif")[Band(1)]
end
# elevationraster = replace(read(elevationraster), elevationraster[10, 10, 1] => -Inf32)

# Population
human_pop = CSV.File(joinpath(workdir, "Population/Population.csv")) |> DataFrame
sugar_cane = CSV.File(joinpath(workdir, "Population/Sugarcane.csv")) |> DataFrame
human_pop.Population .*= 1000
# Plots.plot(human_pop.Year, human_pop.Population)
# Plots.plot(sugar_cane.Year, sugar_cane.Area)

borders = (
    mus=GADM.get("MUS").geom[1],
    reu=GADM.get("REU").geom[1],
)

# Elevation
tiles = getraster(SRTM; bounds=mus_bounds)
dem1 = Raster(tiles[1])
dem2 = Raster(tiles[2])
mus_bounds = (57.1, 57.9), (-20.6, -19.949)
mus_border_selectors = X(Between(mus_bounds[1])), Y(Between(mus_bounds[2])), Band(1)
# Mauritius is right over the split in the tiles
m1 = view(dem1, mus_border_selectors...)
m2 = view(dem2, mus_border_selectors...)
mus_dem = replace_missing(trim(cat(m1, m2; dims=Y); pad=10))
# Plots.plot(mauritius_dem)
reu_bounds = (55.0, 56.0), (-22.0, -20.0)
reu_tile  = getraster(SRTM; bounds=reu_bounds)[1]
reu_border_selectors =  X(Between(reu_bounds[1])), Y(Between(reu_bounds[2])), Band(1)
reu_dem = replace_missing(read(trim(view(Raster(reu_tile), reu_border_selectors...); pad=10)))

# rod_bounds = (63.0, 64.0), (-20.0, -19.0)
# rod_border_selectors =  X(Between(rod_bounds[1])), Y(Between(rod_bounds[2])), Band(1)
# reu_tile  = RasterDataSources.zipurl(SRTM; bounds=rod_bounds)[1]
# reu_tile  = getraster(SRTM; bounds=rod_bounds)[1]
# rod_dem = trim(view(dem3, border_selectors...); pad=10)

dems = (mus=mus_dem, reu=reu_dem)

# Distance
distance_to_coasts = map(dems, borders) do dem, border 
    coast = boolmask(border; to=dem, shape=:line)
    distance_to_coast = nearest_distances(coast)
    mask(distance_to_coast; with=dem)
end

waterways_path = "/home/raf/PhD/Mauritius/Data/water.geojson"
waterways_fc = GeoJSON.read(read(waterways_path))
watermasks = map(dems) do dem 
    boolmask(waterways_fc; to=dem)
end
Plots.plot(watermasks.mus)

distance_to_water = map(dems, watermasks) do dem, watermask
    mask(nearest_distances(watermask); with=dem)
end

Plots.plot(distance_to_water.mus; c=:seaborn_icefire_gradient, size=(1000,1000))
Plots.plot(coasts.reu; c=:seaborn_icefire_gradient, size=(1000,1000))
# Plots.plot!(mus_border; fill=nothing)
# normedelevation = 1 .- elevationraster ./ maximum(elevationraster)
# p1 = distance_to_coast .* distance_to_rivers |> plot;
# p2 = plot(elevationraster)
# p3 = plot(landuse_snapshots[5]; c=:viridis)
# plot(p1, p2, p3; layout=(1, 3))

# Slope
sloperasters = map(d -> slope(d, FD2()), dems)
Plots.plot(sloperasters.reu; clims=(0,0.5))
sloperasters.reu

ps = map(r -> Plots.plot(r; c=:terrain, size=(1000, 1000), clims=(0, 1.0)), sloperasters)
Plots.plot(ps...)

lc_lookup = Dict(lc_categories .=> 0:12)

lc_path = "/home/raf/PhD/Mauritius/Data/Landcover/"
lc_shapes = (
    mus=(
        shapepath=joinpath(lc_path, "mauritius/cla_maurice_fin.shp"),
        crspath=joinpath(lc_path, "mauritius/cla_maurice_fin.prj"),
    ),
    reu=(
        shapepath = joinpath(lc_path, "reunion/cla_run_2014_fin_2975.shp"),
        crspath = joinpath(lc_path, "reunion/cla_run_2014_fin_2975.prj"),
    ),
    # rod=(
        # shapepath = joinpath(lc_path, "rodrigues/cla_rod_fin.shp")
        # crspath = joinpath(lc_path, "rodrigues/cla_rod_fin.prj")
    # )
)

lc_rasters = map(dems, lc_shapes) do dem, (; shapepath, crspath)
    rasterize_lc(dem, shapepath, crspath)
end

# Masks for each land cover
lc_masks = map(lc_rasters) do lc
    rasters = Dict()
    for (k, v) in lc_lookup
        rasters[k] = lc .== v
    end
    rasters
end

Plots.plot(lc_masks.mus["Forest"]; c=:spring)

# plot_lc_makie(lc_rasters.mus)


# Species distributions
#
# GBIF alternatives
# Arctos, neotoma, vertnet

using GBIF, CSV, Plots, DataFrames, Rasters, IntervalSets

species = CSV.File("/home/raf/PhD/Mauritius/mascarine_species.csv") |> DataFrame
for sp in species.Species
    ismissing(sp) && continue
    search = replace(sp, " " => "+")
    run(`chromium https\://www.google.com/search\?q=$(search)`)
    # run(`chromium https\://www.google.com/search\?tbm=isch\&q=$(search)`)
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
