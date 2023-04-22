using Neighborhoods
using Rasters
using Colors
using FileIO
using ImageIO
using GeometryBasics
using GLMakie
using Unitful

using LandscapeChange
using MapRasterization
using GLMakie
includet("common.jl")
includet("map_file_list.jl")
includet("raster_common.jl")
includet("water.jl")

# using JSON3
# using Images
# using GADM
# using GLM
# using ImageCore
# using ColorVectorSpace
# using Dictionaries
# using GeoDataFrames
# using DataFrames, CSV, Tables

# includet("water.jl")
# includet("roads.jl")
# includet("svgs.jl")

# fn = joinpath(datadir, "Selected/Mauritius/Undigitised/mus_landuse_1965_100.png")
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

files = get_map_files()
for island in files
    foreach(file -> choose_categories(first(file); save=false), island)
    yield()
end

# using Plots
img_path = "/home/raf/PhD/Mascarenes/Data/Selected/Mauritius/Undigitised/atlas_18C_land_use_cropped.jpg"
img = load_image(img_path)
# img_path = files.mus[6][1]
# img_path = files.mus.fraser_1835_from_gleadow.filename
output = choose_categories(img_path; save=false, restart=true)
raster_path = splitext(img_path)[1] * ".tif"
r = Raster(raster_path)
ColorSchemes.terrainj
files.mus.atlas_dutch_period.filename


warp_to_raster(files.mus.fraser_1835_from_gleadow.filename, dems.mus;
    object_type=MapRasterization.MapSelection, edit=true,
    guide=(waterways_rivers, waterways_lakes),
)

poly = Polygon([Point2(0.0f0, 1.0f0), Point2(0.0f0, 1.0f0)])
img = load_image(img_path)

using Pkg
Pkg.add(url="https://github.com/rafaqz/MakieDraw.jl")
using MakieDraw
using GLMakie
using ColorSchemes
c = GeometryCanvas{Point2}()
# Click and drag is fast! 

img = rand(1000, 1000)
heatmap!(c.axis, img; colormap=:magma)
# Now its really slow and laggy


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

plotsize = (16, 9) .* 70
plot(RasterStack(slices.mus.timelines.abandonned); size=plotsize, clims=(0, 2))
savefig("abandonned_timeline.png")
plot(RasterStack(map(A -> A .* 2, slices.mus.timelines.cleared)); size=plotsize, clims=(0, 2))
savefig("cleared_timeline.png")
plot(RasterStack(slices.mus.timelines.combined); size=plotsize, clims=(0, 2))
savefig("combined_timeline.png")
plot(slices.mus.files.atlas_1992_agriculture.grouped.urban; c=:viridis, clims=(0, 2))
plot(RasterStack(slices.mus.timelines.urban))

plot(slices.mus.timelines.urban.urban_1992 .& .!(slices.mus.timelines.abandonned.abandonned_1968 .| slices.mus.timelines.cleared.cleared_1968))
plot(slices.mus.timelines.urban.urban_1905)
savefig("urban_1905.png")

map((:mus, :reu), (RasterStack(mus_combined), reu_cleared_timeline), gdal_borders, national_parks) do name, timeline, border, parks
    anim = @animate for A in (values(timeline)..., last(timeline))
        plot(A; legend=false)
        plot!(border; fillalpha=0.5, color=:red)
        plot!(parks; alpha=0, color=:grey, fillalpha=0.5)
        if name == :mus
            for i in 1:3
                class = filter(row -> row.category == i, collect(warped_vegetation))
                c = first(class).color
                color = RGB(c[:r], c[:g], c[:b])
                plot!(class; alpha=0, color, fillalpha=0.5)
            end
        end
    end
    gif(anim, "$(name)_clearing_timeline.gif", fps=0.8)
end
savefig("time_from_ports.png")


using GLMakie
fig = Figure()
ax = Axis(fig[1, 1])
function datashader(limits, pixelarea)
    # return your heatmap data
    # here, I just calculate a sine field as a demo
    xpixels, ypixels = widths(pixelarea)
    xmin, ymin = minimum(limits)
    xmax, ymax = maximum(limits)
    dxmin, dymin = map(max, map(x -> trunc(Int, x), (xmin, ymin)), map(first, axes(A)))
    dxmax, dymax = map(min, map(x -> trunc(Int, x), (xmax, ymax)), map(last, axes(A)))
    # dxpixels = round(Int, (dxmax - dxmin) / (xmax - xmin) * xpixels)
    # dypiyels = round(Int, (dymax - dymin) / (ymax - ymin) * ypixels)
    xs = 
    A[xs, ys]
end
xrange = lift(x -> minimum(x)[1] .. maximum(x)[1], ax.finallimits)
yrange = lift(x -> minimum(x)[2] .. maximum(x)[2], ax.finallimits)
pixels = lift(datashader, ax.finallimits, ax.scene.px_area)
heatmap!(ax, xrange, yrange, pixels,
    xautolimits = false, yautolimits = false)
display(fig)

