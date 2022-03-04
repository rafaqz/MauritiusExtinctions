using CSV
using DataFrames
using GeoJSON
using GADM
# using Plots
using Setfield 
using Shapefile
using Rasters
using Rasters.LookupArrays
using Rasters: set, Between, trim, Band
using GeoInterface
using Plots
using Makie, GLMakie
includet("functions.jl")
years = 1638, 1773, 1835, 1872, 1935, "present"

soilraster = Raster("warpedsoiltypes.tif")[Band(1)]
plot(soilraster)
lakesraster = Raster("warpedlakes.tif")[Band(1)]
elevationraster = Raster("warpedelevation.tif")[Band(1)]
landuse_snapshots = map(1:6) do i
    Raster("warped_landuse_snapshot_$i.tif")[Band(1)]
end
elevationraster = replace(read(elevationraster), elevationraster[10, 10, 1] => -Inf32)


# Population
human_pop = CSV.File("/home/raf/PhD/Mauritius/Data/Population/Population.csv") |> DataFrame
sugar_cane = CSV.File("/home/raf/PhD/Mauritius/Data/Population/Sugarcane.csv") |> DataFrame
human_pop.Population .*= 1000
# Plots.plot(human_pop.Year, human_pop.Population)
# Plots.plot(sugar_cane.Year, sugar_cane.Area)

mus_border = GADM.get("MUS").geom[1]
reu_border = GADM.get("REU").geom[1]
waterways_json = "/home/raf/PhD/Mauritius/Data/osm_rivers.geojson"
waterways = GeoJSON.read(read(waterways_json))
watermask = boolmask(waterways; to=soilraster) .| boolmask(lakesraster)
Plots.plot(watermask)
distance_to_water = mask(nearest_distances(watermask); with=elevationraster[Band(1)])[Band(1)]
Plots.plot(distance_to_water; c=:seaborn_icefire_gradient, size=(1000,1000))
Plots.plot(elevationraster; c=:seaborn_icefire_gradient, size=(1000,1000))
Plots.plot!(watermask; c=:black, legend=:none)
plot!(mus_border; fill=nothing)

i = 1
ps = map(eachindex(years)[2:end]) do i
    year = years[i]
    @show i year
    p = Plots.plot(elevationraster; c=:viridis, legend=:none, ticks=:none, xguide="", yguide="")
    Plots.plot!(p, watermask; c=:blue, legend=:none, opacity=0.5, xguide="", yguide="")
    ss = boolmask(replace(landuse_snapshots[i], 1 => missingval(landuse_snapshots[i])))
    plot!(p, ss; c=:black, legend=:none, opacity=0.5, xguide="", yguide="")
    return p
end

plot(ps...)

savefig("landuse$year.png")

# Elevation
# This is a huge elevation dataset that also happens to be
# in the right projection. Can be downloaded here:
# http://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_DEM/list_5deg.html
# I probably wont actually use this but was for now to make
# the coast and river distance fields.
dem1 = Raster("/home/raf/PhD/Mauritius/DEM/dem_tif_s30e030/s20e055_dem.tif")
dem2 = Raster("/home/raf/PhD/Mauritius/DEM/dem_tif_s30e030/s25e055_dem.tif")
dem3 = Raster("/home/raf/PhD/Mauritius/DEM/dem_tif_s30e060/s20e060_dem.tif")
border_selectors = X(Between(57.1, 57.9)), Y(Between(-20.6, -19.949)), Band(1)
# Mauritius is right over the split in the tiles
m1 = view(dem1, border_selectors...)
m2 = view(dem2, border_selectors...)
mauritius_dem = replace_missing(trim(cat(m1, m2; dims=Y); pad=10))
# Plots.plot(mauritius_dem)
border_selectors =  X(Between(55.0, 56.0)), Y(Between(-20.0, -22.0)), Band(1)
reunion_dem = trim(view(dem2, border_selectors...); pad=10)
border_selectors =  X(63.0..64.0), Y(-20.0..(-19.0)), Band(1)
rodrigues_dem = trim(view(dem3, border_selectors...); pad=10)

# Coast
coast = boolmask(mus_border; to=soilraster, shape=:line)
Plots.plot(coast)
distance_to_coast = nearest_distances(coast)
masked_distance_to_coast = mask(distance_to_coast; with=soilraster[Band(1)])
Plots.plot(masked_distance_to_coast; c=:seaborn_icefire_gradient, size=(1000,1000))
Plots.plot!(mus_border; fill=nothing)
normedelevation = 1 .- elevationraster ./ maximum(elevationraster)
p1 = distance_to_coast .* distance_to_rivers |> plot;
p2 = plot(elevationraster)
p3 = plot(landuse_snapshots[5]; c=:viridis)
plot(p1, p2, p3; layout=(1, 3))

# Slope
sloperaster = slope(mauritius_dem, MaxSlope())
sloperaster = slope(mauritius_dem, FD2())
p1 = Plots.plot(sloperaster; c=:terrain, size=(1000, 1000), clims=(0, 1.0))
savefig("mauritius_slope.png")
p2 = Plots.plot(mauritius_dem; size=(1000, 1000), clims=(0,5))
Plots.plot(p1, p2; size=(2000, 2000))
savefig("mauritius_elevation.png")

lc_categories = [
  "No Data",
  "Continuous urban",
  "Disontinuous urban",
  "Forest",
  "Shrub vegetation",
  "Herbaceaous vegetation",
  "Mangrove",
  "Barren land",
  "Water",
  "Sugarcane",
  "Pasture",
  "",
  "Other cropland",
]
lc_path = "/home/raf/PhD/Mauritius/Data/Landcover/"
mauritius_shapepath = joinpath(lc_path, "mauritius/cla_maurice_fin.shp")
mauritius_crspath = joinpath(lc_path, "mauritius/cla_maurice_fin.prj")
reunion_shapepath = joinpath(lc_path, "reunion/cla_run_2014_fin_2975.shp")
reunion_crspath = joinpath(lc_path, "reunion/cla_run_2014_fin_2975.prj")
rodrigues_shapepath = joinpath(lc_path, "rodrigues/cla_rod_fin.shp")
rodrigues_crspath = joinpath(lc_path, "rodrigues/cla_rod_fin.prj")

mauritius_lc = rasterize_lc(mauritius_dem, mauritius_shapepath, mauritius_crspath);
plot_lc_makie(mauritius_lc)

reunion_lc = rasterize_lc(reunion_dem, reunion_shapepath, reunion_crspath)
plot_lc_makie(reunion_lc)

rodrigues_lc = rasterize_lc(rodrigues_dem, rodrigues_shapepath, rodrigues_crspath)
plot_lc_makie(rodrigues_lc)

x = replace(x -> x == 5 ? x : 0,  mauritius_lc)
plot_lc_makie(x)

lc_shape = Shapefile.Table(mauritius_shapepath)
lc_crs = WellKnownText(readlines(mauritius_crspath)[1])
template = mauritius_dem

lc_df = DataFrame(lc_shape)
lc_raster = Raster(similar(template, Int32); missingval=typemin(Int32))
lc_raster .= typemin(Int32)
lc_raster = read(resample(lc_raster, 100; crs=lc_crs))
fillval = 3
rows = filter(x -> x.ocsol_num == fillval, lc_df)
fillname = first(eachrow(rows)).ocsol_name
# using ProfileView
# @time rasterize!(lc_raster, rows.geometry; fill=fillval)
display(plot_lc(lc_raster))
display(plot_lc_makie(lc_raster))
# hv = lc_raster


# lc_shape = Shapefile.Table(shape_file)
# lc_crs = WellKnownText(readlines(crs_file)[1])
lc_df = DataFrame(lc_shape)
lc_raster = Raster(similar(template, Int32); missingval=typemin(Int32))
lc_raster .= typemin(Int32)
lc_raster = read(resample(lc_raster, 50; crs=lc_crs))
c = Dict(map(=>, lc_categories, 0:12))
# Order of rasterization matters?... (probably should calculate areas?)
fillvals = [
    c["No Data"],
    c["Water"],
    c["Barren land"],
    c["Forest"],
    c["Shrub vegetation"],
    c["Herbaceaous vegetation"],
    c["Mangrove"],
    c["Other cropland"],
    c["Sugarcane"],
    c["Pasture"],
    c["Continuous urban"],
    c["Disontinuous urban"],
]

# @show fillvals
for fillval in fillvals
    rows = filter(x -> x.ocsol_num == fillval, lc_df)
    if length(rows.geometry) > 0
        fillname = first(eachrow(rows)).ocsol_name
        @show fillval, fillname
        rasterize!(lc_raster, rows.geometry; fill=fillval)
    else
        @show fillval 
end end

using Makie, GLMakie
using ProfileView
A = rand(400, 1000)
fig = Figure()
ax, plt = Makie.heatmap(fig[1, 1], A)
display(fig)

filvals = ["string label" for i in 0:10]
    # colormap=cgrad(:cyclic_mygbm_30_95_c78_n256, 10, categorical=true), colorrange=(0, 10),
# )
# ax.aspect = AxisAspect(1)
# Colorbar(fig[1, 2], hm; 
    # ticks=(0:12, lc_categories),
# )


# Arctos, neotoma, vertnet

using GBIF, CSV, Plots, DataFrames, Rasters, IntervalSets
species = CSV.File("/home/raf/PhD/Mauritius/mascarine_species.csv") |> DataFrame
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
RasterDataSources.
prec = Raster(CHELSA{Climate}, :prec; month=1, res)
mascarines_prec = prec[X = Interval(lons...), Y=Interval(lats...)]
Plots.plot(mascarines_prec)
scatter!(points; opacity=0.5, markershape=:circle)
obs = occurrences(, "limit"=>300)

while length(obs) < size(obs)
    @show length(obs) size(obs)
    occurrences!(obs)
end
