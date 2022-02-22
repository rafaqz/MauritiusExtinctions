# Fix the projection of rasters by manually 
# warping them to match another raster
using GLMakie
using Rasters
using Rasters: Band
using Shapefile
using Images
includet("mapfitting.jl")
includet("functions.jl")

# Elevation
# This is a huge elevation dataset that also happens to be
# in the right projection. Can be downloaded here:
# http://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_DEM/list_5deg.html
# I probably wont actually use this but was for now to make
# the coast and river distance fields.
# dem1 = Raster("/home/raf/PhD/Mauritius/DEM/dem_tif_s30e030/s20e055_dem.tif")
# dem2 = Raster("/home/raf/PhD/Mauritius/DEM/dem_tif_s30e030/s25e055_dem.tif")
# border_selectors =  X(Between(57.1, 57.9)), Y(Between(-20.6, -19.949)), Band(1)
# m1 = view(dem1, border_selectors...)
# m2 = view(dem2, border_selectors...)
# dem = replace_missing(trim(cat(m1, m2; dims=Y); pad=5))
# # Plots.plot(dem)

# # Land use files
# elevpath = "/home/raf/PhD/Mauritius/Data/Norder/LS factor/DEM/DEM100x100_Resample.img"
# lakespath = "/home/raf/PhD/Mauritius/Data/Norder/LS factor/Lakes/lakes_all.shp"
# soiltypespath = "/home/raf/PhD/Mauritius/Data/Norder/K factor/SoilK.shp"
# rainfallpath = "/home/raf/PhD/Mauritius/Data/Norder/R factor/r_annual.img"
# landusedir = "/home/raf/PhD/Mauritius/Data/Norder/C factor/"

# rainfallraster = Raster(rainfallpath)[Band(1)]
# # Elevation is a slightly larger raster for some reason
# # `crop` doesn't work because the index is slightly different
# elevationraster = Raster(elevpath; missingval=-3.4028235f38)[Band(1), X(1:520)]
# elevationraster = Raster(elevpath; missingval=-3.4028235f38)[Band(1)]
# # Plots.plot(elevation)

# lakes_shapes = Shapefile.Handle(lakespath)
# lakesraster = zeros(Union{Int32,Missing}, dims(rainfallraster))
# lakesraster .= missing
# for i in eachindex(lakes_shapes.shapes)[1:end-2]
#     rasterize!(lakesraster, lakes_shapes.shapes[i]; fill=i)
# end

# soiltypes_shapes = Shapefile.Handle(soiltypespath)
# soilraster = copy(lakesraster) .= missing
# for i in eachindex(soiltypes_shapes.shapes)
#     rasterize!(soilraster, soiltypes_shapes.shapes[i]; fill=i)
# end
# elevationraster = mask(elevationraster; with=soilraster)
# # Plots.plot(soilraster)
# Plots.plot(elevationraster)
# # Plots.plot(lakesraster)

# year = 1935
# years = 1638, 1773, 1835, 1872, 1935, "present"
# landuse_shapefiles = map(years) do year
#     path = joinpath(landusedir, string(year, ".shp"))
#     Shapefile.Handle(path)
# end
# landuse_snapshots = map(landuse_shapefiles, years) do shapefile, year
#     landuse = copy(lakesraster) .= missing
#     # The forested/cleared order swaps after the first three files
#     shapes = year in (1773, 1835) ? shapefile.shapes : reverse(shapefile.shapes)
#     for (n, shape) in enumerate(shapes)
#         rasterize!(landuse, shape; fill=n)
#     end
#     return landuse
# end

# landuse_snapshots
# # Plots.plot(Plots.plot.(landuse_snapshots; clims=(0, 2), c=:viridis)...; size=(2000,2000))

# i = 1
# ps = map(eachindex(years)[2:end]) do i
#     year = years[i]
#     @show i year
#     p = Plots.plot(elevationraster; c=:viridis, legend=:none, ticks=:none, xguide="", yguide="")
#     ss = boolmask(replace(landuse_snapshots[i], 1 => missingval(landuse_snapshots[i])))
#     plot!(p, ss; c=:black, legend=:none, opacity=0.5, xguide="", yguide="")
#     return p
# end

# wsoilraster, wlakesraster, welevationraster, wlandusesnapshots... = 
#     MakieRasters.manualwarp(soilraster, lakesraster, elevationraster, landuse_snapshots...; to=dem)
# write("warpedsoiltypes.tif", wsoilraster)
# write("warpedlakes.tif", wlakesraster)
# write("warpedelevation.tif", welevationraster)
# for i in 1:length(landuse_snapshots)
#     write("warped_landuse_snapshot_$i.tif", wlandusesnapshots[i])
# end


img = load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page33_mauritius_vegetation.png")
img = load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page36_reunion_vegetation.png")
img = load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page145_reunion_phases.png")
img = load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page145_reunion_phases_colored.png")
img = load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page157_mauritius_settlements.png")
img = load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page159_mauritius_remnants.png")
img = load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page166_mauritius_fodies.png")
img = load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page176_reunion_settlements.png")
img = load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page183_reunion_remnants.png")
img = load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page194_rodrigues_settlements.png")
img = load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page252_mauritius_kestrel.png")
img = load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page265_reunion_bulbul.png")

mauritius_remnants_categories = (0.7058823529411765, 0.2823529411764706)

remnants = Float64.((x -> x.val).(Gray.(img)))
rast = Raster(rotr90(remnants), (X, Y))
Makie.heatmap(parent(rast))

points = MakieRasters.manualinput(rast)
points = MakieRasters.manualinput(rast; points)

counts = Dict{Float64,Int}()
foreach(rast) do x
    ismissing(x) && return
    if haskey(counts, x)
        counts[x] += 1
    else
        counts[x] = 1
    end
end
counts
for (x, c) in counts
    (c < 6000 || x == 0) && delete!(counts, x)
end
categories = Tuple(keys(counts))
categories = (0.23137254901960785, 0.8, 0.45098039215686275, 0.9411764705882353, 0.615686274509804)

warped = MakieRasters.manualwarp(rast; to=reunion_dem)
warped1 = replace_missing(warped, 1.0)

cleaned = clean_categories(warped1; categories, neighborhood=Moore{2}(), missingval=0.0)
cleaned = clean_categories(cleaned; categories, neighborhood=Moore{2}(), missingval=0.0)

cleaned = reverse(cleaned; dims=2)
cleaned = replace(cleaned, NaN => 1.0)
replaced = replace(x -> x < 0.1 ? 0.0 : 1.0, warped1)
p = Plots.plot(replace_missing(reunion_dem) ./ maximum(reunion_dem))
p = Plots.plot!(p, replace_missing(cleaned); opacity=0.5)


heatmap(parent(mask(cleaned; with=dem)))
warpeb = rebuild(reverse(warped; dims=2); dims=dims(dem))

dims(warped)
dims(dem)
