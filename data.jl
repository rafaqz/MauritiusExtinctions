using CSV
using DataFrames
using GeoJSON
# using GADM
# using Plots
using Setfield 
using Shapefile
using Rasters
using Rasters.LookupArrays
using Rasters: set, Between, trim
using Plots

include("functions.jl")

# Land use files
elevpath = "/home/raf/PhD/Mauritius/Data/Norder/LS factor/DEM/DEM100x100_Resample.img"
lakespath = "/home/raf/PhD/Mauritius/Data/Norder/LS factor/Lakes/lakes_all.shp"
soiltypespath = "/home/raf/PhD/Mauritius/Data/Norder/K factor/SoilK.shp"
rainfallpath = "/home/raf/PhD/Mauritius/Data/Norder/R factor/r_annual.img"
landusedir = "/home/raf/PhD/Mauritius/Data/Norder/C factor/"
waterways_json = "/home/raf/PhD/Mauritius/Data/osm_rivers.geojson"

rainfallraster = Raster(rainfallpath)[Band(1)]
# Elevation is a slightly larger raster for some reason
# `crop` doesn't work because the index is slightly different
elevationraster = Raster(elevpath; missingval=-3.4028235f38)[Band(1), X(1:520)]
# Plots.plot(elevation)

lakes_shapes = Shapefile.Handle(lakespath)
lakesraster = Raster{Int}(undef, dims(rainfall); missingval=0)
for i in eachindex(lakes_shapes.shapes)[1:end-1]
    rasterize!(lakesraster, lakes_shapes.shapes[i]; fill=i)
end


# soiltypes_shapes = Shapefile.Handle(soiltypespath)
soilraster = lakesraster .* 0
for i in eachindex(soiltypes_shapes.shapes)
    rasterize!(soilraster, soiltypes_shapes.shapes[i]; fill=i)
end
elevationraster = mask(elevationraster; with=soilraster)

year = 1935
years = 1638, 1773, 1835, 1872, 1935, "present"
landuse_shapefiles = map(years) do year
    path = joinpath(landusedir, string(year, ".shp"))
    Shapefile.Handle(path)
end

landuse_snapshots = map(landuse_shapefiles, years) do shapefile, year
    landuse = zeros(Int, dims(elevation))
    # The forested/cleared order swaps after the first three files
    shapes = year in (1773, 1835) ? shapefile.shapes : reverse(shapefile.shapes)
    for (n, shape) in enumerate(shapes)
        rasterize!(landuse, shape; fill=n)
    end
    return landuse
end
# Plots.plot(plot.(landuse_snapshots; clims=(0, 2), c=:viridis)...; size=(2000,2000))
# warpedsoilraster = manual_warp(soilraster; to=dem)

# wsoilraster, wlakesraster, welevationraster, wlandusesnapshots... = 
    # manualwarp(soilraster, lakesraster, elevationraster, landuse_snapshots...; to=dem)
# TODO: fix Rasters.jl write so this isn't needed
# _asint32(A) = Int32.(replace_missing(A, typemin(Int32)))
# write("warpedsoiltypes.tif", _asint32(wsoilraster))
# write("warpedlakes.tif", _asint32(wlakesraster))
# write("warpedelevation.tif", welevationraster)
# for i in 1:length(landuse_snapshots)
#     write("landuse_snapshot_$i.tif", _asint32(landuse_snapshots[i]))
# end
# As = Makie.heatmap(parent(wlakesraster))

# display(Makie.heatmap(map(parent, dims(wlandusesnapshots[1]))..., parent(wlandusesnapshots[1])))

soilraster = Raster("warpedsoiltypes.tif")
lakesraster = Raster("warpedlakes.tif")
elevationraster = Raster("warpedelevation.tif")
Plots.plot(soilraster)
landuse_snapshots = map(1:6) do i
    Raster("landuse_snapshot_$i.tif")
end
Plots.plot(landuse_snapshots[4])

# Population
human_pop = CSV.File("/home/raf/PhD/Mauritius/Data/Population/Population.csv") |> DataFrame
sugar_cane = CSV.File("/home/raf/PhD/Mauritius/Data/Population/Sugarcane.csv") |> DataFrame
human_pop.Population .*= 1000
# Plots.plot(human_pop.Year, human_pop.Population)
# Plots.plot(sugar_cane.Year, sugar_cane.Area)

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
Plots.plot(dem)

# Plots.plot(dem, size=(1000,1000))
# Plots.plot(dem; c=:gist_earth, size=(1000,1000))
# Plots.plot(dem; c=:terrain, size=(1000,1000))
# Plots.plot(dem; c=:cubehelix, size=(1000,1000))
# Plots.plot(dem; c=:batlow, size=(1000,1000))
# Plots.plot(dem; c=:seaborn_icefire_gradient, size=(1000,1000))

# Rivers
# waterways = GeoJSON.read(read(waterways_json))
# waterways.crs
# rivers = boolmask(waterways; to=dem)
# Plots.plot!(rivers; c=:seaborn_icefire_gradient)
# Plots.plot!(waterways; c=:green)
# distance_to_rivers = mask(nearest_distances(rivers); with=dem)
# Plots.plot(distance_to_rivers; c=:gist_earth, size=(1000,1000))
# Plots.plot(distance_to_rivers; c=:batlow, size=(1000,1000))
# Plots.plot(distance_to_rivers; c=:seaborn_icefire_gradient, size=(1000,1000))
# Plots.plot!(waterways; c=:green)
# plot!(mauritius_border; fill=nothing)
# plot!(rivers; c=:seaborn_icefire_gradient)

# Coast
# mauritius_border = GADM.get("MUS").geom[1]
# Plots.plot!(mauritius_border; fill=nothing)
# coast = boolmask(mauritius_border; to=dem, shape=:line)
# Plots.plot(coast)
# distance_to_coast = nearest_distances(coast)
# masked_distance_to_coast = mask(distance_to_coast; with=dem)
# Plots.plot(masked_distance_to_coast; c=:seaborn_icefire_gradient, size=(1000,1000))
# Plots.plot!(waterways; c=:green)
# Plots.plot!(mauritius_border; fill=nothing)
# distance_to_coast .* distance_to_rivers .* dem |> plot

# Slope
# sloperaster = slope(elevation, MaxSlope())
# sloperaster = slope(elevation, FD2())
# p1 = Plots.plot(sloperaster; c=:terrain, size=(1000, 1000), clims=(0, 1.0))
# savefig("mauritius_slope.png")
# p2 = Plots.plot(dem; size=(1000, 1000), clims=(0,5))
# Plots.plot(p1, p2; size=(2000, 2000))
# savefig("mauritius_elevation.png")
