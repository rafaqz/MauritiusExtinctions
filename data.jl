using CSV
using GLM
using DataFrames
using GADM
using GeoJSON
using GeoFormatTypes
using GeoInterface
using Plots
using Rasters
using Setfield 
using Shapefile
using NearestNeighbors
using Rasters.LookupArrays
using Rasters: set, Between, trim
using Plots


# Land use files
elevpath = "/home/raf/PhD/Mauritius/Data/Norder/LS factor/DEM/DEM100x100_Resample.img"
lakespath = "/home/raf/PhD/Mauritius/Data/Norder/LS factor/Lakes/lakes_all.shp"
soiltypespath = "/home/raf/PhD/Mauritius/Data/Norder/K factor/SoilK.shp"
rainfallpath = "/home/raf/PhD/Mauritius/Data/Norder/R factor/r_annual.img"
landusedir = "/home/raf/PhD/Mauritius/Data/Norder/C factor/"
waterways_json = "/home/raf/PhD/Mauritius/Data/osm_rivers.geojson"

elevation = Raster(elevpath; missingval=-3.4028235f38)[Band(1)]
rainfall = Raster(rainfallpath)

lakes_shapes = Shapefile.Handle(lakes_path)
lakesraster = Raster{Int}(undef, dims(r); missingval=0)
for i in eachindex(lakes_shapes.shapes)[1:end-1]
    rasterize!(lakesraster, lakes.shapes[i]; fill=i)
end

soiltypes_shapes = Shapefile.Handle(soiltypespath)
soilraster = lakesraster .* 0
for i in eachindex(soiltypes_shapes.shapes)
    rasterize!(soilraster, soiltypes.shapes[i]; fill=i)
end


years = 1638, 1773, 1835, 1872, 1935, "present"
landuse_shapefiles = map(years) do year
    path = joinpath(landuse_dir, string(year, ".shp"))
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
plot(plot.(landuse_snapshots; clims=(0, 2), c=:viridis)...; size=(2000,2000))

# DEM
mauritius_border = GADM.get("MUS").geom[1]

human_pop = CSV.File("/home/raf/PhD/Mauritius/Data/Population/Population.csv") |> DataFrame
sugar_cane = CSV.File("/home/raf/PhD/Mauritius/Data/Population/Sugarcane.csv") |> DataFrame
human_pop.Population .*= 1000
plot(human_pop.Year, human_pop.Population)
plot(sugar_cane.Year, sugar_cane.Area)

dem1 = Raster("/home/raf/PhD/Mauritius/Data/DEM/dem_tif_s30e030/s20e055_dem.tif")
dem2 = Raster("/home/raf/PhD/Mauritius/Data/DEM/dem_tif_s30e030/s25e055_dem.tif")
border_selectors =  X(Between(57.1, 57.9)), Y(Between(-20.6, -19.949)), Band(1)
m1 = view(dem1, border_selectors...)
m2 = view(dem2, border_selectors...)
dem = replace_missing(trim(cat(m1, m2; dims=Y); pad=5))
plot(dem, size=(1000,1000))
plot(dem; c=:gist_earth, size=(1000,1000))
plot(dem; c=:terrain, size=(1000,1000))
plot(dem; c=:cubehelix, size=(1000,1000))
plot(dem; c=:batlow, size=(1000,1000))
plot(dem; c=:seaborn_icefire_gradient, size=(1000,1000))

# Rivers
waterways = GeoJSON.read(read(waterways_json))
waterways.crs

rivers = boolmask(waterways; to=dem)
# plot!(rivers; c=:seaborn_icefire_gradient)
plot!(waterways; c=:green)
plot!(mauritius_border; fill=nothing)

distance_to_rivers = mask(nearest_distances(rivers); with=dem)
plot(distance_to_rivers; c=:gist_earth, size=(1000,1000))
plot(distance_to_rivers; c=:batlow, size=(1000,1000))
plot(distance_to_rivers; c=:seaborn_icefire_gradient, size=(1000,1000))
plot!(waterways; c=:green)
plot!(mauritius_border; fill=nothing)
# plot!(rivers; c=:seaborn_icefire_gradient)

coast = boolmask(mauritius_border; to=dem, shape=:polygon)
plot(coast)
distance_to_coast = mask(nearest_distances(coast); with=dem)
plot(distance_to_coast; c=:seaborn_icefire_gradient, size=(1000,1000))
plot!(waterways; c=:green)
plot!(mauritius_border; fill=nothing)
distance_to_coast .* distance_to_rivers .* dem |> plot
# plot!(rivers; c=:seaborn_icefire_gradient)

# using ProfileView
# @profview slope_ = slope(dem, MaxSlope())
p1 = plot(slope_; c=:terrain, size=(1000, 1000), clims=(0, 1.0))
savefig("mauritius_slope.png")
p2 = plot(dem; size=(1000, 1000), clims=(0,5))
plot(p1, p2; size=(2000, 2000))
savefig("mauritius_elevation.png")
