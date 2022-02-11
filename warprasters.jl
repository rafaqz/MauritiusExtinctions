# Fix the projection of rasters by manually 
# warping them to match another raster

include("mapfitting.jl")

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

# Land use files
elevpath = "/home/raf/PhD/Mauritius/Data/Norder/LS factor/DEM/DEM100x100_Resample.img"
lakespath = "/home/raf/PhD/Mauritius/Data/Norder/LS factor/Lakes/lakes_all.shp"
soiltypespath = "/home/raf/PhD/Mauritius/Data/Norder/K factor/SoilK.shp"
rainfallpath = "/home/raf/PhD/Mauritius/Data/Norder/R factor/r_annual.img"
landusedir = "/home/raf/PhD/Mauritius/Data/Norder/C factor/"

rainfallraster = Raster(rainfallpath)[Band(1)]
# Elevation is a slightly larger raster for some reason
# `crop` doesn't work because the index is slightly different
elevationraster = Raster(elevpath; missingval=-3.4028235f38)[Band(1), X(1:520)]
# Plots.plot(elevation)

lakes_shapes = Shapefile.Handle(lakespath)
lakesraster = zeros(Union{Int32,Missing}, dims(rainfallraster))
lakesraster .= missing
for i in eachindex(lakes_shapes.shapes)[1:end-2]
    rasterize!(lakesraster, lakes_shapes.shapes[i]; fill=i)
end

soiltypes_shapes = Shapefile.Handle(soiltypespath)
soilraster = copy(lakesraster) .= missing
for i in eachindex(soiltypes_shapes.shapes)
    rasterize!(soilraster, soiltypes_shapes.shapes[i]; fill=i)
end
elevationraster = mask(elevationraster; with=soilraster)
Plots.plot(soilraster)
Plots.plot(elevationraster)
Plots.plot(lakesraster)

year = 1935
years = 1638, 1773, 1835, 1872, 1935, "present"
landuse_shapefiles = map(years) do year
    path = joinpath(landusedir, string(year, ".shp"))
    Shapefile.Handle(path)
end
landuse_snapshots = map(landuse_shapefiles, years) do shapefile, year
    landuse = copy(lakesraster) .= missing
    # The forested/cleared order swaps after the first three files
    shapes = year in (1773, 1835) ? shapefile.shapes : reverse(shapefile.shapes)
    for (n, shape) in enumerate(shapes)
        rasterize!(landuse, shape; fill=n)
    end
    return landuse
end
Plots.plot(Plots.plot.(landuse_snapshots; clims=(0, 2), c=:viridis)...; size=(2000,2000))

wsoilraster, wlakesraster, welevationraster, wlandusesnapshots... = 
    manualwarp(soilraster, lakesraster, elevationraster, landuse_snapshots...; to=dem)
write("warpedsoiltypes.tif", wsoilraster)
write("warpedlakes.tif", wlakesraster)
write("warpedelevation.tif", welevationraster)
for i in 1:length(landuse_snapshots)
    write("warped_landuse_snapshot_$i.tif", wlandusesnapshots[i])
end
