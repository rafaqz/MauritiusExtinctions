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
rainfallpath = "/home/raf/PhD/Mauritius/Data/Norder/R factor/r_annual.img"

rainfallraster = Raster(rainfallpath)[Band(1)]
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
# landusedir = "/home/raf/PhD/Mauritius/Data/Norder/C factor/"
# landuse_shapefiles = map(years) do year
#     path = joinpath(landusedir, string(year, ".shp"))
#     Shapefile.Handle(path)
# end
# landuse_snapshots = map(landuse_shapefiles, years) do shapefile, year
#     landuse = zeros(Union{Int32,Missing}, dims(rainfallraster)) .= missing
#     # The forested/cleared order swaps after the first three files
#     shapes = year in (1773, 1835) ? shapefile.shapes : reverse(shapefile.shapes)
#     for (n, shape) in enumerate(shapes)
#         rasterize!(landuse, shape; fill=n)
#     end
#     return landuse
# end
# Plots.plot(landuse_snapshots[1])
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
        veg=load(joinpath(datadir, "LostLand/Maps/page33_mauritius_vegetation_colored.png")),
        phase=load(joinpath(datadir, "LostLand/Maps/page132_mauritius_phases_colored.png")),
        rem=load(joinpath(datadir, "LostLand/Maps/page159_mauritius_remnants_colored.png")),
        # kestrel=load(joinpath(datadir, "LostLand/Maps/page252_mauritius_kestrel_colored.png")),
    )),
    reu=map(rotr90, (
        veg=load(joinpath(datadir, "LostLand/Maps/page36_reunion_vegetation_colored.png")),
        phase=load(joinpath(datadir, "LostLand/Maps/page145_reunion_phases_colored.png")),
        rem=load(joinpath(datadir, "LostLand/Maps/page183_reunion_remnants_colored.png")),
        # bulbul=load(joinpath(datadir, "LostLand/Maps/page265_reunion_bulbul_colored.png")),
    ))
)

island_image_classes = (
    mus=(
        veg=(
            red = "semi-dry evergreen forest",
            green = "open dry palm-rich woodland",
            blue = "wet forest",
            cyan = "Pandanus swamp",
            magenta = "mossy rainforest",
            yellow = "mangrove",
            dark_red = "wetland vegetation",
        ),
        phase=(
            magenta = "deforested before 1807",
            dark_red = "deforested 1807-1835",
            blue = "deforested 1835-1910",
            dark_blue = "deforested 1910-1947",
            green = "deforested 1947-1970",
            red = "deforested since 1970",
            yellow = "scrub with native remnants",
            dark_green = "native vegetation",
            cyan = "swamp",
        ),
        rem=(
            red = "at least 20% native plants",
            green = "50%+ native plants",
            blue = "reservoirs",
            cyan = "cleared",
        ),
    ),
    reu=(
        veg=(
            dark_cyan = "Palm-bejoin savanna",
            magenta = "Palm-rich dry forest",
            red = "Lowland mixed evergreen semi-dry tropical forest",
            yellow = "Lowland mixed evergreen tropical rainforest", 
            cyan = "Transitional mixed evergreen tropical rainforest",
            green = "Montane mixed evergreen subtropical rainforest",
            dark_blue = "Tree-heather formations",
            dark_green = "Acacia-dominated montane rainforest",
            dark_red = "Very wet screw-pine formations",
            blue = "Upper montane temperate heaths",
            dark_yellow = "Wetland vegetation",
            dark_magenta = "Unvegetated recent lava flows",
        ),
        phase=(
            blue = "17th century",
            red = "18th century",
            green = "19th century",
            yellow = "20th century",
            cyan = "not cleared/colonised",
        ),
        rem=(
            red = "Lowland mixed evergreen semi-dry forest",
            green = "Lowland mixed evergreen tropical forest",
            blue = "Transitional mixed evergreen tropical rainforest",
            yellow = "Montane mixed evergreen subtropical rainforest",
            magenta = "Tree-heather formations",
            cyan = "Acacia-dominated montane rainforest",
            dark_red = "Very wet screw-pine formations",
            dark_green = "Upper montane temperate heaths",
            dark_blue = "Wetland vegetation",
            dark_yellow = "Largely unvegetated recent lava flows",
            dark_magenta = "Cleared",
        ),
        # bulbul=load(joinpath(datadir, "LostLand/Maps/page265_reunion_bulbul_colored.png")),
    )
)

# m_stl = load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page157_mauritius_settlements_colored.png") |> rotr90
# m_fodies = load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page166_mauritius_fodies_colored.png") |> rotr90
# re_stl = load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page176_reunion_settlements.png") |> rotr90
# ro_stl = load("/home/raf/PhD/Mauritius/Data/LostLand/Maps/page194_rodrigues_settlements.png") |> rotr90

_color(name::String) = RGB{N0f8}((Colors.color_names[name] ./ 255)...)
colors_color_names = ("red1", "green1", "blue1", "cyan", "magenta", "yellow")
main_color_names = Symbol.(strip.(isnumeric, colors_color_names))
dark_color_names = Symbol.(string.(Ref("dark_"), main_color_names))
main_colors = NamedTuple{main_color_names}(RGBA{N0f8}.(_color.(colors_color_names)))
dark_colors = NamedTuple{dark_color_names}(RGBA{N0f8}.(RGB.(values(main_colors)) .* 0.498)) # GIMP rounds round, Colors.jl rounds up...
color_lookup = merge(main_colors, dark_colors)

island_image_colors = map(island_image_classes) do image_classes
    map(image_classes) do classes
        colors = map(c -> color_lookup[c], keys(classes))
    end
end

templates = map(dems, borders) do dem, border
    mask(dem; with=boolmask(border; to=dem, shape=:polygon, boundary=:touches))
    dem
end
# Plots.plot(templates.mus)

# Mass point selection for all category filtered images
# foreach(keys(island_images), island_images) do i, images
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
# island_points = map(keys(island_images), island_images) do i, images
#     map(keys(images)) do k
#         path = joinpath(outputdir, "WarpPoints/$(i)_$(k)_warp_points.csv")
#         isfile(path) ? CSV.File(path) : missing
#     end |> NamedTuple{keys(images)}
# end |> NamedTuple{keys(island_images)}

# # Run the category cleaning filter
# cleaned_island_images = map(island_images, island_image_colors) do images, image_colors
#     map(images, image_colors) do A, colors
#         clean_categories(A; 
#             categories=colors, neighborhood=Moore{3}(), missingval=RGBA{N0f8}(1.0, 1.0, 1.0, 1.0)
#         )
#     end
# end
# # Apply warping from raw images to cleaned images
# island_image_rasters = map(keys(island_images), island_images, cleaned_island_images, island_points) do i, images, cimages, layerpoints
#     map(keys(images), images, cimages, layerpoints) do k, A1, A2, points
#         rs = RasterUtils.applywarp(A1, A2; 
#             template=templates[i], points=DataFrame(points), missingval=RGBA{N0f8}(1.0,1.0,1.0)
#         )
#         (raw=rs[1], cleaned=rs[2])
#     end |> NamedTuple{keys(images)}
# end |> NamedTuple{keys(island_images)}
# # Re-classify cleaned color categories as integers
# island_rasters = map(island_image_rasters, island_image_colors) do ir, ic
#     map(ir, ic) do rs, colors
#         A = rs.cleaned
#         classify(A, map(=>, colors, UInt8.(eachindex(colors)))...; others=UInt8(0), missingval=UInt8(0))
#     end
# end
# # Write images and plots to disk
# foreach(keys(island_rasters), island_rasters, island_image_classes, dems, borders) do i, ir, ic, dem, border
#     foreach(keys(ir), ir, ic) do k, A, classes
#         grad = cgrad([map(c -> color_lookup[c], keys(classes))...])
#         p = Plots.plot(A; c=grad, clims=(1, length(grad)))
#         Plots.plot!(p, border; fill=nothing)
#         savefig(p, joinpath(outputdir, "$(i)_$(k)_cleaned.png"))
#         p = Plots.plot(dem ./ maximum(skipmissing(dem)))
#         Plots.plot!(p, A ./ maximum(A); c=:spring, opacity=0.6)
#         Plots.plot!(p, border; fill=nothing)
#         savefig(p, joinpath(outputdir, "$(i)_$(k)_cleaned_dem.png"))
#         write(joinpath(outputdir, "$(i)_$(k)_cleaned.tif"), A) 
#     end
# end

# Reload rasters from disk
island_rasters = map(keys(island_images), island_images) do i, ii
    map(keys(ii)) do k
        read(Raster(joinpath(outputdir, "$(i)_$(k)_cleaned.tif")))
    end |> NamedTuple{keys(ii)}
end |> NamedTuple{keys(island_images)}

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

# Plots.plot(island_rasters.mus.veg)
# keys(island_rasters.mus)

# map(1:6)

# Vegetation phases as total cleared for each time step
veg_phases = map(island_rasters, (mus=0:6, reu=0:4)) do island, phaseids
    map(phaseids) do id
        classify(island.phase, 1..id => true; others=false, missingval=missing)
    end
end

cleared_rasters = map(island_rasters) do i
    cleared = maximum(i.rem)
    (x -> x == cleared).(replace_missing(i.rem))
end

cleared_rem_plots = map(veg_phases, cleared_rasters) do vp, rem
    Plots.plot(
        Plots.plot((!).(vp[end]); c=:spring, title="Not cleared"), 
        Plots.plot((!).(rem); c=:spring, title="Remnant")
    )
end
cleared_rem_plots.reu
p = cleared_rem_plots.mus
savefig(p, "not_cleared_mus.png")


