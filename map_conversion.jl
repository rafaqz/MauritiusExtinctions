using Stencils
using Rasters
using Colors
using FileIO
using ImageIO
using GeometryBasics
using GLMakie
using Unitful
import GeoInterface as GI

using MapRasterization
using GLMakie
includet("common.jl")
includet("map_file_list.jl")
includet("raster_common.jl")
includet("water.jl")

swap_ext(path, newext) = string(splitext(path)[1], newext)

# using JSON3
# using Images
# using GADM
# using GLM
# using ImageCore
# using ColorVectorSpace
# using Dictionaries
# using GeoDataFrames
# using DataFrames, CSV, Tables

# includet("roads.jl")
# includet("svgs.jl")

# output = choose_categories(fn; save=true, restart=true)

#=
L. Maillard map drawn 1845-52
Categories:
Cultives tropicales
Forets
Laves Modernes
Laves anciennes et Saverunes?a??
Jardinage
=#

files = define_map_files()
for island in files
    foreach(file -> choose_categories(first(file); save=false), island)
    yield()
end
#

img_path = "/home/raf/PhD/Mascarenes/maps/Rodrigues/Gade_1985_native.png"
run(`gimp $img_path`)
img = load_image(img_path)
output = MapRasterization.CategorySelector(img; ncategories=5)
output.figure
summary = MapRasterization.CategorySummary(output);
ms = MapRasterization.MapSelection(output)
JSON3.write(swap_ext(img_path, ".json"), ms)
# run(`gimp $img_path`)
# using Plots
img_path = "/home/raf/PhD/Mascarenes/Data/Selected/Mauritius/Undigitised/atlas_19C_land_use.jpg"
# run(`gimp $img_path`)
img_path = joinpath(datadir, "Selected/Mauritius/Undigitised/mus_landuse_1965_100_hi_c.png")
img = load_image(img_path)
# size(img)
# img2 = img[1:32:2000, 1:32:2000]
# size(img2)
# using SegmentAnything
# x = SegmentAnything.generate(SegmentAnything.MaskGenerator(), img2)
# y = PythonCall.pyconvert(Vector{Dict}, x)
# segs = map(d -> d["segmentation"], y) |> stack
# heatmap(sum(segs; dims=3)[:, :, 1])
Makie.set_theme!(theme_black())
# features = GeoJSON.read(read("mus_landuse_1965_100_polygons.json"))
output = MapRasterization.CategorySelector(img; polygons=features, ncategories=11)
summary = MapRasterization.CategorySummary(output);
ms = MapRasterization.MapSelection(output)
JSON3.write(swap_ext(img_path, ".json"), ms)
GeoJSON.write("mus_landuse_1965_100_polygons.json", summary.polygons)
# rasterize(count, summary.polygons, size=(100, 100))
# display(output)

# warpselector = MapRasterization.WarpSelector(img; point_table, template=dems.mus, guide=borders.mus)
# CSV.write("mus_landuse_1965_100_warp_points.csv", warpselector.warper.point_table[])
point_table = CSV.File("mus_landuse_1965_100_warp_points.csv")
multipolygons = GI.convert.(MultiPolygon, (GI.geometry(f) for f in features))
warper = MapRasterization.Warper(; point_table, template=dems.mus, poly=1)
warped = warper(multipolygons; missingval=zero(RGBAf))

map(getproperty(:name), GI.getfeature(features))
lc_names = map(getproperty(:name), GI.getfeature(features))

order = [
    "Sugar",
    "Scrub",
    "Forrest_Natural",
    "Forrest_plantation",
    "Vegetables",
    "Savannah",
    "Swamps",
    "Reservoirs",
    "Builtp",
    "Tea",
    "Rock",
]
rast = fill!(similar(dems.mus, Union{Int,Missing}), missing)
for s in order
    n = findfirst(==(s), lc_names)
    rasterize!(rast, warped[n]; fill=n)
end
mask!(rast; with=dems.mus)
raster_path = swap_ext(img_path, ".tif")
write(raster_path, rast)
Rasters.rplot(rast; colormap=:batlow)

# img1 = img[map(x -> 1:2:x, size(img))...]
# output = MapRasterization.CategorySelector(img; ncategories=11)
# GeoJSON.write("mus_landuse_1965_100_points.json", summary.points)
# polygons = GeoJSON.read(read("mus_landuse_1965_100_polygons.json"))
# GeoInterface.trait(polygons)

# using SegmentAnything
# SegmentAnything.unsafe_empty_cache()

# warp_to_raster(img_path, dems.mus;
#     # object_type=MapRasterization.MapSelection, 
#     object_type=MapRasterization.CategoryShapes,
#     edit=true,
#     guide=(waterways_rivers, waterways_lakes),
# )


# Click and drag is fast! 

# choose_categories(files.reu.cadet_invasives.filename)
# open_output(files.mus["atlas_19C_land_use"]).settings.category_name

# Readable json
# io = IOBuffer()
# pretty_json_path = splitext(img_path)[1] * "_pretty.json"
# JSON3.pretty(io, output)
# write(pretty_json_path, String(take!(io)))

# Warping

# Warp
# map(get_map_files(), map(replace_missing, dems), borders, road_lines) do island, template, border, roads
#     map(island) do (filepath, poly)
#         warp_to_raster(filepath, template .* 0; poly, edit=false, guide=(roads, border, waterways_rivers, waterways_lakes))
#     end
# end
slices = make_raster_slices(masks, lc_categories)

using Tyler
using Extents, TileProviders
Tyler.Map(Extent(X=(1, 10), Y=(0, 10)); provider=Google())
web = resample(slices.mus.timelines.lc[end-1:end]; crs=EPSG(3857))
Makie.heatmap!(web[2]; transparency=true, colormap=(:red, 0.2))

plotsize = (16, 9) .* 70
using Plots
Plots.plot(web[1])
Plots.plot(slices.mus.timelines.urban ; size=plotsize, clims=(0, 2))
savefig("abandonned_timeline.png")
Plots.plot(map(A -> A .* 2, slices.mus.timelines.cleared); size=plotsize, clims=(0, 2))
savefig("cleared_timeline.png")
Plots.plot(slices.mus.timelines.lc[end-2:end]; size=plotsize, clims=(0, 5))
savefig("combined_timeline.png")
Plots.plot(slices.mus.files.atlas_1992_agriculture.grouped.urban; c=:viridis, clims=(0, 2))
Plots.plot(slices.mus.timelines.urban)



plot(slices.mus.timelines.urban.urban_1992 .& .!(slices.mus.timelines.abandonned.abandonned_1968 .| slices.mus.timelines.cleared.cleared_1968))
plot(slices.mus.timelines.urban.urban_1905)
savefig("urban_1905.png")
timeline = slices.mus.timelines.lc
name = :mus
border = gdal_borders.mus

# map((:mus, :reu), (mus_combined, reu_cleared_timeline), gdal_borders) do name, timeline, border
    anim = Plots.@animate for A in (values(timeline)..., last(timeline))
        Plots.plot(A; legend=false, clims=(0, 5), c=:batlow)
        # Plots.plot!(border; fillalpha=0.0)
        # Plots.plot!(parks; alpha=0, color=:grey, fillalpha=0.5)
        if name == :mus
            for i in 1:3
                class = filter(row -> row.category == i, collect(mus_native_veg_poly))
                c = first(class).color
                color = RGB(c["r"], c["g"], c["b"])
                Plots.plot!(GI.geometry(class); alpha=0, color, fillalpha=0.5)
            end
        end
    end
    Plots.gif(anim, "$(name)_clearing_timeline.gif", fps=0.8)
# end
savefig("time_from_ports.png")


using Tyler
using MakieDraw
using Extents
using TileProviders
using GLMakie
using GeometryBasics
using GeoJSON
using JSON3
includet("raster_common.jl")

geoms = GeoJSON.read(read("reunion_1790.json"))
figure = Figure()
axis = Axis(figure[1:10, 1:10])
tyler = Tyler.Map(Extents.extent(dems.reu); provider=Google(), figure, axis, scale=2)
canvas = MakieDraw.GeometryCanvas{Polygon}(canvas; #[Polygon([Point(1.0, 2.0), Point(1.0, 2.0)])]; 
    figure=tyler.figure,
    axis=tyler.axis,
    # properties=(;landcover=Int[]),
    propertynames=(:landcover,),
    poly_kw=(; colorrange=(1, 2), transparency=true),
)
display(tyler.figure)
GeoJSON.write("reunion_1790.json", canvas)
