using CSV
using DataFrames
using GeoJSON
using GADM
# using Plots
using Setfield 
using Shapefile
using Rasters
using Rasters.LookupArrays
using Rasters: set, Between, trim
using Plots

include("functions.jl")
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

mauritius_border = GADM.get("MUS").geom[1]
waterways_json = "/home/raf/PhD/Mauritius/Data/osm_rivers.geojson"
waterways = GeoJSON.read(read(waterways_json))
watermask = boolmask(waterways; to=soilraster) .| boolmask(lakesraster)
plot(watermask)
distance_to_water = mask(nearest_distances(watermask); with=elevationraster[Band(1)])[Band(1)]
Plots.plot(distance_to_water; c=:seaborn_icefire_gradient, size=(1000,1000))
Plots.plot(elevationraster; c=:seaborn_icefire_gradient, size=(1000,1000))
Plots.plot!(watermask; c=:black, legend=:none)
plot!(mauritius_border; fill=nothing)

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
border_selectors =  X(Between(57.1, 57.9)), Y(Between(-20.6, -19.949)), Band(1)
m1 = view(dem1, border_selectors...)
m2 = view(dem2, border_selectors...)
dem = replace_missing(trim(cat(m1, m2; dims=Y); pad=5))
# Plots.plot(dem)

# Coast
coast = boolmask(mauritius_border; to=soilraster, shape=:line)
Plots.plot(coast)
distance_to_coast = nearest_distances(coast)
masked_distance_to_coast = mask(distance_to_coast; with=soilraster[Band(1)])
Plots.plot(masked_distance_to_coast; c=:seaborn_icefire_gradient, size=(1000,1000))
Plots.plot!(mauritius_border; fill=nothing)
normedelevation = 1 .- elevationraster ./ maximum(elevationraster)
p1 = distance_to_coast .* distance_to_rivers |> plot;
p2 = plot(elevationraster)
p3 = plot(landuse_snapshots[5]; c=:viridis)
plot(p1, p2, p3; layout=(1, 3))

# Slope
sloperaster = slope(elevation, MaxSlope())
sloperaster = slope(elevation, FD2())
p1 = Plots.plot(sloperaster; c=:terrain, size=(1000, 1000), clims=(0, 1.0))
savefig("mauritius_slope.png")
p2 = Plots.plot(dem; size=(1000, 1000), clims=(0,5))
Plots.plot(p1, p2; size=(2000, 2000))
savefig("mauritius_elevation.png")


# Elevation
# This is a huge elevation dataset that also happens to be
# in the right projection. Can be downloaded here:
# http://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_DEM/list_5deg.html
# I probably wont actually use this but was for now to make
# the coast and river distance fields.
dem1 = Raster("/home/raf/PhD/Mauritius/DEM/dem_tif_s30e030/s20e055_dem.tif")
dem2 = Raster("/home/raf/PhD/Mauritius/DEM/dem_tif_s30e030/s25e055_dem.tif")
border_selectors =  X(Between(55.0, 56.0)), Y(Between(-20.0, -22.0)), Band(1)
reunion_dem = trim(view(dem2, border_selectors...))
plot(reunion_dem)
