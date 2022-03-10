# Fix the projection of rasters by manually 
# warping them to match another raster
using GLMakie
using Rasters
using Rasters: Band
using Shapefile
using DataFrames
using Images
using CSV
using Colors, ColorVectorSpace
includet("mapfitting.jl")
includet("functions.jl")

# Land use files
# elevpath = "/home/raf/PhD/Mauritius/Data/Norder/LS factor/DEM/DEM100x100_Resample.img"
# lakespath = "/home/raf/PhD/Mauritius/Data/Norder/LS factor/Lakes/lakes_all.shp"
# soiltypespath = "/home/raf/PhD/Mauritius/Data/Norder/K factor/SoilK.shp"
# rainfallpath = "/home/raf/PhD/Mauritius/Data/Norder/R factor/r_annual.img"
# landusedir = "/home/raf/PhD/Mauritius/Data/Norder/C factor/"

# rainfallraster = Raster(rainfallpath)[Band(1)]
# Elevation is a slightly larger raster for some reason
# `crop` doesn't work because the index is slightly different
# elevationraster = Raster(elevpath; missingval=-3.4028235f38)[Band(1), X(1:520)]
# Plots.plot(rainfallraster)

# lakes_shapes = Shapefile.Handle(lakespath)
# lakesraster = zeros(Union{Int32,Missing}, dims(rainfallraster))
# lakesraster .= missing
# for i in eachindex(lakes_shapes.shapes)[1:end-2]
#     rasterize!(lakesraster, lakes_shapes.shapes[i]; fill=i)
# end

# soiltypes_shapes = Shapefile.Handle(soiltypespath)
# soilraster = replace_missing(similar(rainfallraster, Int)) .= missing
# for i in eachindex(soiltypes_shapes.shapes)
#     rasterize!(soilraster, soiltypes_shapes.shapes[i]; fill=i)
# end

# elevationraster = mask(elevationraster; with=soilraster)
# Plots.plot(soilraster)
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
# Plots.plot(Plots.plot.(landuse_snapshots; clims=(0, 2), c=:viridis)...; size=(2000,2000))

# i = 1
# ps = map(eachindex(years)[2:end]) do i
#     year = years[i]
#     @show i year
#     p = Plots.plot(elevationraster; c=:viridis, legend=:none, ticks=:none, xguide="", yguide="")
#     ss = boolmask(replace(landuse_snapshots[i], 1 => missingval(landuse_snapshots[i])))
#     Plots.plot!(p, ss; c=:black, legend=:none, opacity=0.5, xguide="", yguide="")
#     return p
# end
# Plots.plot(ps...)

# wsoilraster, wlakesraster, welevationraster, wrainfalraster, wlandusesnapshots... = 
#     MakieRasters.manualwarp(soilraster, lakesraster, elevationraster, rainfallraster, landuse_snapshots...; to=dem)
# write("warpedsoiltypes.tif", wsoilraster)
# write("warpedlakes.tif", wlakesraster)
# write("warpedrainfall.tif", wlakesraster)
# write("warpedelevation.tif", welevationraster)
# for i in 1:length(landuse_snapshots)
#     write("warped_landuse_snapshot_$i.tif", wlandusesnapshots[i])
# end

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

templates = map(dems, borders) do dem, border
    mus=boolmask(border; to=dem, shape=:line)
end

# Mass point selection for all category filtered images
foreach(keys(island_images), island_images) do i, images
    foreach(keys(images), images) do k, A
        path = joinpath(outputdir, "WarpPoints/$(i)_$(k)_warp_points.csv")
        points = isfile(path) ? CSV.File(path) : nothing
        warp_points = RasterUtils.select_common_points(Float64.(Gray.(A));
            template=templates[i], missingval=missing, points=points
        )
        CSV.write(joinpath(outputdir, "WarpPoints/$(i)_$(k)_warp_points.csv"), warp_points)
    end
end

# Load previously selected points
island_points = map(keys(island_images), island_images) do i, images
    map(keys(images)) do k
        path = joinpath(outputdir, "WarpPoints/$(i)_$(k)_warp_points.csv")
        @show path
        isfile(path) ? CSV.File(path) : missing
    end |> NamedTuple{keys(images)}
end |> NamedTuple{keys(island_images)}

# Run the category cleaning filter
cleaned_island_images = map(island_images) do images
    map(images) do A
        clean_categories(A; 
            categories=basic_colors, neighborhood=Moore{2}(), missingval=RGBA{N0f8}(1.0, 1.0, 1.0, 1.0)
        )
    end
end
Makie.image(cleaned_island_images.reu.rem)

# Apply warping to all images
island_rasters = map(keys(island_images), island_images, cleaned_island_images, island_points) do i, images, cimages, points
    map(keys(images), images, cimages) do k, A1, A2
        rs = RasterUtils.manualwarp(A, A2; 
            template=templates[i], points, missingval=RGBA{N0f8}(1.0,1.0,1.0)
        )
        (raw=rs[1], cleaned=rs[2])
    end |> NamedTuple{keys(images)}
end |> NamedTuple{keys(island_images)}

# Macaque distribution
mus_macaques = load(joinpath(datadir, "Macaques/macaque_distribution.png")) |> rotr90
mus_macaques_mask = load(joinpath(datadir, "Macaques/macaque_distribution_mask.png")) |> rotr90
mus_macaques_bool = Float64.(Gray.(mus_macaques)) .< 0.9
mus_macaque_warp_points = RasterUtils.select_common_points(m_macaques; 
    template, missingval=RGB{N0f8}(1.0,1.0,1.0),
    points=DataFrame(mus_macaque_warp_points)
)
CSV.write(joinpath(outputdir, "WarpPoints/macaque_warp_points.csv"), mus_macaque_warp_points)
mus_macaque_warp_points = CSV.File(joinpath(datadir, "Macaques/macaque_warp_points.csv"))
w_mus_macaques, w_mus_macaques_bool = RasterUtils.manualwarp(mus_macaques, mus_macaques_bool; 
    template, points=mus_macaque_warp_points, missingval=RGB{N0f8}(1.0,1.0,1.0)
)
Plots.plot(Float64.(Gray.(w_mus_macaques)))
Plots.plot(Float64.(Gray.(mus_macaques_bool)); c=:spring)
Plots.plot!(borders.mus; fill=nothing)
write(joinpath(datadir, "Macaques/macaques_mask.tif"), wm_macaques_bool)

# Kestrel distribution
# Has multiple files so we need to do it separately
kestrel_warp_points = CSV.File(joinpath(outputdir("WarpPoints/mus_kestrel_warp_points.csv"))
w_mus_kestrel_2001, w_mus_kestrel_cleaned, w_mus_kestrel = RasterUtils.manualwarp(
    mus_kestrel_2001_mask, mus_kestrel_cleaned, mus_kestrel; 
    template=templates.mus, points=kestrel_warp_points, missingval=RGB{N0f8}(1.0,1.0,1.0)
)
Plots.plot(Float64.(Gray.(w_mus_kestrel)), c=:spring)
Plots.plot!(Float64.(Gray.(w_mus_kestrel_2001)); opacity=0.5, c=:spring)
Plots.plot!(borders.mus; fill=nothing)

# kestrel_1830_1860 = m_kestrel_cleaned .==  
# kestrel_1930_1950 = m_kestrel_cleaned .== 
# kestrel_1970_1980 = m_kestrel_cleaned .== 

write("/home/raf/PhD/Mauritius/Data/Macaques/kestrel.tif", w_m_kestrel_cleaned)
