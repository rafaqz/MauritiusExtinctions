using Neighborhoods
using MapRasterization
using JSON3
using Images
using GLMakie
using GADM
using Colors
using ImageIO
using FileIO
using GLM
using ImageCore
using ColorVectorSpace
using Dictionaries
using GeometryBasics
using Rasters
using GeoDataFrames

# json = "/home/raf/PhD/Mascarenes/Data/Selected/Mauritius/Undigitised/1835_fraser_composite_etsy.json"
# using Mmap
# JSON3.read(Mmap.mmap(json))

includet("svgs.jl")
includet("map_file_list.jl")

#=
L. Maillard map drawn 1845-52
Categories:
Cultives tropicales
Forets
Laves Modernes
Laves anciennes et Saverunes???
Jardinage
=#

files = get_map_files()

for island in files
    foreach(file -> choose_categories(first(file); save=false), island)
    yield()
end

# choose_categories(files.mus["atlas_1992_agriculture"].filename)
# open_output(files.mus["atlas_19C_land_use"]).settings.category_name

# Readable json
# io = IOBuffer()
# pretty_json_path = splitext(img_path)[1] * "_pretty.json"
# JSON3.pretty(io, output)
# write(pretty_json_path, String(take!(io)))

# Warping
using DataFrames, CSV, Tables
includet("raster_common.jl")
includet("water.jl")
includet("roads.jl")

# Warp
map(get_map_files(), map(replace_missing, dems), borders, road_lines) do island, template, border, roads
    map(island) do (filepath, poly)
        warp_to_raster(filepath, template .* 0; poly, edit=false, guide=(roads, border, waterways_rivers, waterways_lakes))
    end
end

slices = make_raster_slices(masks)

plot(RasterStack(slices.reu.timelines.cleared))
plot(RasterStack(slices.mus.timelines.abandonned))
plot(RasterStack(slices.mus.timelines.cleared))
plot(RasterStack(slices.mus.timelines.urban))
plot(RasterStack(slices.mus.timelines.combined); c=:viridis, clims=(0, 2))
plot(slices.mus.files.atlas_1992_agriculture.grouped.urban; c=:viridis, clims=(0, 2))

plot(slices.mus.timelines.urban.urban_1992 .& .!(slices.mus.timelines.abandonned.abandonned_1968 .| slices.mus.timelines.cleared.cleared_1968))
plot(slices.mus.timelines.urban.urban_1905)
savefig("urban_1905.png")

map((:mus, :reu), (RasterStack(mus_combined), reu_cleared_timeline), gdal_borders, national_parks) do name, timeline, border, parks
    anim = @animate for A in values(timeline)
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
