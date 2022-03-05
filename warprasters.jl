# Fix the projection of rasters by manually 
# warping them to match another raster
using GLMakie
using Rasters
using Rasters: Band
using Shapefile
using Images
using CSV
using Colors, ColorVectorSpace
includet("mapfitting.jl")
includet("functions.jl")

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
Plots.plot(rainfallraster)

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
# elevationraster = mask(elevationraster; with=soilraster)
# Plots.plot(soilraster)
Plots.plot(elevationraster)
# Plots.plot(lakesraster)

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

landuse_snapshots
# Plots.plot(Plots.plot.(landuse_snapshots; clims=(0, 2), c=:viridis)...; size=(2000,2000))

i = 1
ps = map(eachindex(years)[2:end]) do i
    year = years[i]
    @show i year
    p = Plots.plot(elevationraster; c=:viridis, legend=:none, ticks=:none, xguide="", yguide="")
    ss = boolmask(replace(landuse_snapshots[i], 1 => missingval(landuse_snapshots[i])))
    Plots.plot!(p, ss; c=:black, legend=:none, opacity=0.5, xguide="", yguide="")
    return p
end
Plots.plot(ps...)

wsoilraster, wlakesraster, welevationraster, wrainfalraster, wlandusesnapshots... = 
    MakieRasters.manualwarp(soilraster, lakesraster, elevationraster, rainfallraster, landuse_snapshots...; to=dem)
write("warpedsoiltypes.tif", wsoilraster)
write("warpedlakes.tif", wlakesraster)
write("warpedrainfall.tif", wlakesraster)
write("warpedelevation.tif", welevationraster)
for i in 1:length(landuse_snapshots)
    write("warped_landuse_snapshot_$i.tif", wlandusesnapshots[i])
end

# Cleaning re-colored maps
island_images = (
    mus=map(rotr90, (
        veg=load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page33_mauritius_vegetation_colored.png"),
        kestrel=load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page252_mauritius_kestrel_colored.png"),
        rem=load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page159_mauritius_remnants_colored.png"),
    )),
    reu=map(rotr90, (
        veg=load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page36_reunion_vegetation_colored.png"),
        phase=load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page145_reunion_phases_colored.png"),
        rem=load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page183_reunion_remnants_colored.png"),
        bulbul=load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page265_reunion_bulbul_colored.png"),
    ))
)

# m_stl = load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page157_mauritius_settlements_colored.png") |> rotr90
# m_fodies = load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page166_mauritius_fodies_colored.png") |> rotr90
# re_stl = load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page176_reunion_settlements.png") |> rotr90
# ro_stl = load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page194_rodrigues_settlements.png") |> rotr90

color(name::String) = RGB{N0f8}((Colors.color_names[name] ./ 255)...)

colors = color.(("white", "red", "blue", "green", "yellow", "cyan", "magenta"))
white, red, blue, green, yellow, cyan, magenta = maincolors = RGBA{N0f8}.(colors)
_, darkred, darkblue, darkgreen, darkyellow, darkcyan, darkmagenta = halfcolors = RGBA{N0f8}.(colors ./ 2)

basic_colors = (maincolors..., halfcolors...)

templates = (
    mus=boolmask(mus_border; to=mauritius_dem, shape=:line),
    reu=boolmask(reu_border; to=reunion_dem, shape=:line),
)

cleaned_island_images = map(island_images) do images
    map(images) do A
        clean_categories(A; 
            categories=basic_colors, neighborhood=Moore{2}(), missingval=RGBA{N0f8}(1.0, 1.0, 1.0, 1.0)
        )
    end
end
Makie.image(cleaned_island_images.reu.rem)

# Mass point selection for all images
foreach(keys(island_images), island_images) do i, images
    foreach(keys(images), images) do k, A
        warp_points = RasterUtils.select_common_points(Float64.(Gray.(A)); template=templates[i])
        CSV.write("/home/raf/PhD/Mauritius/Data/LostLand/$(i)_$(k)_warp_points.csv", warp_points)
    end
end

# Mass warping for all images
island_rasters = map(keys(island_images), island_images, cleaned_island_images) do i, images, cimages
    map(keys(images), images, cimages) do k, A1, A2
        warp_points = CSV.File("/home/raf/PhD/Mauritius/Data/LostLand/$(i)_$(k)_warp_points.csv")
        rs = RasterUtils.manualwarp(A, A2; 
            template=templates[i], points=warp_points, 
            missingval=RGBA{N0f8}(1.0,1.0,1.0)
        )
        (raw=rs[1], cleaned=rs[2])
    end |> NamedTuple{keys(images)}
end |> NamedTuple{keys(cleaned_island_images)}

# Macaque distribution
m_macaques = load("/home/raf/PhD/Mauritius/Data/Macaques/macaque_distribution.png") |> rotr90
m_macaques_mask = load("/home/raf/PhD/Mauritius/Data/Macaques/macaque_distribution_mask.png") |> rotr90
m_macaques_bool = Float64.(Gray.(m_macaques)) .< 0.9
# macaque_warp_points = RasterUtils.select_common_points(m_macaques_reference; template)
# CSV.write("/home/raf/PhD/Mauritius/Data/Macaques/macaque_warp_points.csv", macaque_warp_points)
macaque_warp_points = CSV.File("/home/raf/PhD/Mauritius/Data/Macaques/macaque_warp_points.csv")
w_m_macaques, w_m_macaques_bool = RasterUtils.manualwarp(m_macaques, m_macaques_bool; 
    template, points=macaque_warp_points, missingval=RGB{N0f8}(1.0,1.0,1.0)
)
Plots.plot(Float64.(Gray.(w_m_macaques)))
Plots.plot(Float64.(Gray.(m_macaques_bool)); c=:spring)
Plots.plot!(mauritius_border; fill=nothing)
write("/home/raf/PhD/Mauritius/Data/Macaques/macaques_mask.tif", wm_macaques_bool)

# Kestrel distribution
m_kestrel = load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page252_mauritius_kestrel_colored.png") |> rotr90
m_kestrel_2001 = load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page252_mauritius_kestrel_2001.png") |> rotr90
m_kestrel_2001_mask = Float64.(Gray.(m_kestrel_2001)) .< 0.5

# CSV.write("/home/raf/PhD/Mauritius/Data/LostLand/kestrel_warp_points.csv", kestrel_warp_points)
kestrel_warp_points = CSV.File("/home/raf/PhD/Mauritius/Data/LostLand/kestrel_warp_points.csv")
w_m_kestrel_2001, w_m_kestrel_cleaned, w_m_kestrel = RasterUtils.manualwarp(m_kestrel_2001_mask, m_kestrel_cleaned, m_kestrel; 
    template, points=kestrel_warp_points, missingval=RGB{N0f8}(1.0,1.0,1.0)
)
Plots.plot(Float64.(Gray.(w_m_kestrel)), c=:spring)
Plots.plot!(Float64.(Gray.(w_m_kestrel_2001)); opacity=0.5, c=:spring)
Plots.plot!(mauritius_border; fill=nothing)

# kestrel_1830_1860 = m_kestrel_cleaned .==  
# kestrel_1930_1950 = m_kestrel_cleaned .== 
# kestrel_1970_1980 = m_kestrel_cleaned .== 

write("/home/raf/PhD/Mauritius/Data/Macaques/kestrel.tif", w_m_kestrel_cleaned)
