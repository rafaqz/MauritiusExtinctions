# Fix the projection of rasters by manually 
# warping them to match another raster
using Rasters
using Rasters.LookupArrays
using Rasters: set, Between, trim, Band
using RasterUtils
using Shapefile
using DataFrames
using Images
using CSV
using GeoInterface
using Colors, ColorVectorSpace
using Statistics
using GLMakie

includet("raster_common.jl")
includet("functions.jl")
includet("lost_land_images.jl")

templates = map(dems, borders) do dem, border
    mask(dem; with=boolmask(border; to=dem, shape=:polygon, boundary=:touches))
    dem
end

# Norder land use files
elevpath = "$datadir/Norder/LS factor/DEM/DEM100x100_Resample.img"
lakespath = "$datadir/Norder/LS factor/Lakes/lakes_all.shp"
soiltypespath = "$datadir/Norder/K factor/SoilK.shp"
rainfallpath = "$datadir/Norder/R factor/r_annual.img"

raw_rainfallraster = Raster(rainfallpath, crs=EPSG(3337))[Band(1)]
# Elevation is a slightly larger raster for some reason
# `crop` doesn't work because the index is slightly different
raw_elevationraster = Raster(elevpath; missingval=-3.4028235f38, crs=EPSG(3337))[Band(1), X(1:520)]

lakes_shapes = Shapefile.Handle(lakespath)
raw_lakesraster = zeros(Int32, dims(raw_rainfallraster); missingval=typemin(Int32))
raw_lakesraster .= typemin(Int32)
for i in eachindex(lakes_shapes.shapes)[1:end-2]
    rasterize!(raw_lakesraster, lakes_shapes.shapes[i]; fill=i)
end
soiltypes_shapes = Shapefile.Handle(soiltypespath)
raw_soilraster = replace_missing(similar(raw_rainfallraster, Int32), typemin(Int32)) .= typemin(Int32)
for i in eachindex(soiltypes_shapes.shapes)
    rasterize!(raw_soilraster, soiltypes_shapes.shapes[i]; fill=i)
end
plot(raw_soilraster)
raw_elevationraster = mask(raw_elevationraster; with=raw_soilraster)

landusedir = "$datadir/Norder/C factor/"
landuse_shapefiles = map(years) do year
    path = joinpath(landusedir, string(year, ".shp"))
    Shapefile.Handle(path)
end
raw_landuse_rasters = map(landuse_shapefiles, years) do shapefile, year
    landuse = zeros(Int32, dims(raw_rainfallraster); missingval=typemin(Int32)) .= typemin(Int32)
    # The forested/cleared order swaps after the first three files
    shapes = year in (1773, 1835) ? shapefile.shapes : reverse(shapefile.shapes)
    for (n, shape) in enumerate(shapes)
        rasterize!(landuse, shape; fill=n)
    end
    return landuse
end
raw_landuse_stack = RasterStack(raw_landuse_rasters; keys=lc_year_keys)

raw_norder_stack = merge(RasterStack((
    soiltypes=raw_soilraster,
    lakes=raw_lakesraster,
    elevation=raw_elevationraster,
    rainfall=raw_rainfallraster,
)), raw_landuse_stack)

warped_norder_stack = resample(raw_norder_stack; to=dems.mus, filename="test.tif")
norder_dir = mkpath(joinpath(outputdir, "Norder/"))
write(string(norder_dir, "/"), warped_norder_stack; ext=".tif")
plot(RasterStack(norder_dir)[Band(1)]; c=:viridis)

m_stl = load("$datadir/LostLand/Maps/page157_mauritius_settlements_colored.png") |> rotr90
m_fodies = load("$datadir/LostLand/Maps/page166_mauritius_fodies_colored.png") |> rotr90
re_stl = load("$datadir/LostLand/Maps/page176_reunion_settlements.png") |> rotr90
ro_stl = load("$datadir/LostLand/Maps/page194_rodrigues_settlements.png") |> rotr90
# Plots.plot(templates.mus)

m_stl = load("$datadir/LostLand/Maps/page157_mauritius_settlements_colored.png") |> rotr90

# Lost Land map images from pdf

_color(name::String) = RGB{N0f8}((Colors.color_names[name] ./ 255)...)
colors_color_names = ("red1", "green1", "blue1", "cyan", "magenta", "yellow")
main_color_names = Symbol.(strip.(isnumeric, colors_color_names))
dark_color_names = Symbol.(string.(Ref("dark_"), main_color_names))
main_colors = NamedTuple{main_color_names}(RGBA{N0f8}.(_color.(colors_color_names)))
dark_colors = NamedTuple{dark_color_names}(RGBA{N0f8}.(RGB.(values(main_colors)) .* 0.498)) # GIMP rounds round, Colors.jl rounds up...
color_lookup = merge(main_colors, dark_colors)

lostland_images = map(lostland_image_paths) do island
    map(rotr90 ∘ load, island)
end
lostland_image_colors = map(lostland_image_classes) do image_classes
    map(image_classes) do classes
        colors = map(c -> color_lookup[c], keys(classes))
    end
end

# Mass point selection for all category filtered images
# foreach(keys(lostland_images), lostland_images) do i, images
#     foreach(keys(images), images) do k, A
#         path = joinpath(outputdir, "WarpPoints/$(i)_$(k)_warp_points.csv")
#         points = isfile(path) ? CSV.File(path) : nothing
#         warp_points = RasterUtils.select_common_points(Float64.(Gray.(A));
#             template=templates[i], missingval=missing, points=points
#         )
#         CSV.write(joinpath(outputdir, "WarpPoints/$(i)_$(k)_warp_points.csv"), warp_points)
#     end
# end

# Load previously selected points
lostland_points = map(namedkeys(lostland_images), lostland_images) do i, images
    map(namedkeys(images)) do k
        path = joinpath(outputdir, "WarpPoints/$(i)_$(k)_warp_points.csv")
        isfile(path) ? CSV.File(path) : missing
    end
end

# Run the category cleaning filter
cleaned_lostland_images = map(lostland_images, lostland_image_colors) do images, image_colors
    map(images, image_colors) do A, colors
        clean_categories(A; 
            categories=colors, neighborhood=Moore{3}(), missingval=RGBA{N0f8}(1.0, 1.0, 1.0, 1.0)
        )
    end
end

# Apply warping from raw images to cleaned images
lostland_image_rasters = map(keys(lostland_images), lostland_images, cleaned_lostland_images, lostland_points) do i, images, cimages, layerpoints
    map(keys(images), images, cimages, layerpoints) do k, A1, A2, points
        rs = RasterUtils.applywarp(A1, A2; 
            template=templates[i], points=DataFrame(points), missingval=RGBA{N0f8}(1.0,1.0,1.0)
        )
        (raw=rs[1], cleaned=rs[2])
    end |> NamedTuple{keys(images)}
end |> NamedTuple{keys(lostland_images)}
# Re-classify cleaned color categories as integers
lostland_rasters = map(lostland_image_rasters, lostland_image_colors) do ir, ic
    map(ir, ic) do rs, colors
        A = rs.cleaned
        classify(A, map(=>, colors, UInt8.(eachindex(colors)))...; others=UInt8(0), missingval=UInt8(0))
    end
end
# Write images and plots to disk

foreach(keys(lostland_rasters), lostland_rasters, lostland_image_classes, dems, borders) do i, ir, ic, dem, border
    foreach(keys(ir), ir, ic) do k, A, classes
        dir = joinpath(outputdir, "LostLand", string(i))
        plotdir = joinpath(outputdir, "LostLand/Plots", string(i))
        mkpath(dir)
        mkpath(plotdir)
        grad = cgrad([map(c -> color_lookup[c], keys(classes))...])
        p = Plots.plot(A; c=grad, clims=(1, length(grad)))
        Plots.plot!(p, border; fill=nothing)
        savefig(p, joinpath(outputdir, "LostLand/Plots", string(i), "$(k).png"))
        p = Plots.plot(dem ./ maximum(skipmissing(dem)))
        Plots.plot!(p, A ./ maximum(A); c=:spring, opacity=0.6)
        Plots.plot!(p, border; fill=nothing)
        savefig(p, joinpath(outputdir, "LostLand/Plots", string(i), "$(k)_dem.png"))
        write(joinpath(outputdir, "LostLand", string(i), "$(k).tif"), A) 
    end
end

# Macaque distribution
# mus_macaques = load(joinpath(datadir, "Macaques/macaque_distribution.png")) |> rotr90
# mus_macaques_mask = load(joinpath(datadir, "Macaques/macaque_distribution_mask.png")) |> rotr90
# mus_macaques_bool = Float64.(Gray.(mus_macaques)) .< 0.9
# mus_macaque_warp_points = RasterUtils.select_common_points(m_macaques; 
#     template, missingval=RGB{N0f8}(1.0,1.0,1.0),
#     points=DataFrame(mus_macaque_warp_points)
# )
# CSV.write(joinpath(outputdir, "WarpPoints/macaque_warp_points.csv"), mus_macaque_warp_points)
# mus_macaque_warp_points = CSV.File(joinpath(datadir, "Macaques/macaque_warp_points.csv"))
# w_mus_macaques, w_mus_macaques_bool = RasterUtils.manualwarp(mus_macaques, mus_macaques_bool; 
#     template=templates.mus,
#     # points=DataFrame(mus_macaque_warp_points),
#     missingval=RGB{N0f8}(1.0,1.0,1.0)
# )
# Plots.plot(Float64.(Gray.(w_mus_macaques)))
# Plots.plot(Float64.(Gray.(mus_macaques_bool)); c=:spring)
# Plots.plot!(borders.mus; fill=nothing)
# write(joinpath(datadir, "Macaques/macaques_mask.tif"), wm_macaques_bool)

# Kestrel distribution
# Has multiple files so we need to do it separately
# kestrel_warp_points = CSV.File(joinpath(outputdir("WarpPoints/mus_kestrel_warp_points.csv")))
# w_mus_kestrel_2001, w_mus_kestrel_cleaned, w_mus_kestrel = RasterUtils.manualwarp(
#     mus_kestrel_2001_mask, mus_kestrel_cleaned, mus_kestrel; 
#     template=templates.mus, points=kestrel_warp_points, missingval=RGB{N0f8}(1.0,1.0,1.0)
# )
# Plots.plot(Float64.(Gray.(w_mus_kestrel)), c=:spring)
# Plots.plot!(Float64.(Gray.(w_mus_kestrel_2001)); opacity=0.5, c=:spring)
# Plots.plot!(borders.mus; fill=nothing)

# Plots.plot(lostland_rasters.mus.veg)
# keys(lostland_rasters.mus)

# map(1:6)

rem = load("$datadir/LostLand/Maps/page159_mauritius_remnants.png") |> rotr90
RasterUtils.selectcolors(rem)
