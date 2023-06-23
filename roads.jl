
includet("raster_common.jl")

# This is not one file, but a series of "mus_roads_19XX"
roads_paths = (;
    mus = joinpath(datadir, "Generated/Roads/mus/mus_roads.tif"),
    reu = joinpath(datadir, "Generated/Roads/reu/reu_roads.tif"),
    rod = joinpath(datadir, "Generated/Roads/reu/rod_roads.tif"),
)

# Empty roads for rodrigues
# rodroads = falses(axes(dems.rod))
# write(roads_paths.rod, rodroads)

road_series = map(roads_paths, dems[keys(roads_paths)]) do path, dem
    series = map(A -> Bool.(A), RasterSeries(path, Ti))
    mask(series; with=dem, missingval=missing)
end

# Plots.plot(road_series.mus)
# Plots.plot(road_series.reu)


# using Dates, DataStructures
# using DimensionalData.LookupArrays
# using GeometryBasics
# using GeoInterface
# using GeoInterfaceRecipes
# using GeoJSON
# GeoInterfaceRecipes.@enable_geo_plots GeoJSON.Feature

# includet("raster_common.jl")
# includet("waynames.jl")

# To roads
# function namelike(names, feature)
#     haskey(feature.properties, "name") || return false
#     return any(occursin(feature.properties["name"]), names)
# end
# function roadin(feature, roadlist)
#     hasproperty(feature, :highway) || return false
#     feature.highway in ("trunk", "primary", "secondary", "tertiary", "unclassified", "track", "residential") || return false
#     return hasproperty(feature, :ref) && getproperty(feature, :ref) in roadlist ||
#         hasproperty(feature, :name) && getproperty(feature, :name) in roadlist ||
#         hasproperty(feature, Symbol("@id")) && getproperty(feature, Symbol("@id")) in roadlist
# end
# select_ways(ways, way_names) =
#     GeoInterface.FeatureCollection([f for f in ways if roadin(f, way_names)])
# function plot_ways(ways, way_names, dem, title="")
#     selected_ways = select_ways(ways, way_names)
#     plot(dem ./ maximum(skipmissing(dem)), legend=false)
#     # plot((dems.mus ./ maximum(skipmissing(dems.mus)))[X=57.4..58.6, Y=(-20.35)..(-20.15)])
#     # plot!(watermasks.mus; color=:blue, legend=false)
#     # plot!(deforestation.mus[:by_1835]; c=:viridis, opacity=0.4)a
#     # plot!(selected_ways; color=:red, title, legend=false)
#     plot!(selected_ways; title, legend=false)
# end
# function rasterize_ways(ways, way_names; to)
#     selected_ways = select_ways(ways, way_names)
#     rasterize(selected_ways; to, fill=true)
# end

# @time ways = map(island_keys) do ik
#     json_path = joinpath(datadir, "Roads", "$(ik)_ways.geojson")
#     GeoJSON.read(Base.read(json_path))
# end

# way_names_sequences = (;
#     mus = (
#         1667=>mus_way_names_1667,
#         1725=>mus_way_names_1725,
#         1764=>mus_way_names_1764,
#         1813=>mus_way_names_1813,
#         1880=>mus_way_names_1880,
#     ),
#     reu = (
#         1752=>reu_way_names_1752,
#         1804=>reu_way_names_1804,
#     )
# )

# selected_ways = select_ways(ways.mus, way_names_sequences.mus[4][end])

# @time road_lines = map(ways) do island_ways
#     lines = GeometryBasics.LineString[]
#     for (i, way) in enumerate(island_ways)
#         println(i)
#         geom = GeoInterface.geometry(way)
#         if geom isa GeoJSON.LineString
#             push!(lines, GeoInterface.convert(GeometryBasics.LineString, geom))
#         end
#     end
#     lines
# end
# using GLMakie
# GLMakie.lines(road_lines.mus)

# way_masks = map(ways, way_names_sequences, dems) do w, seq, dem
#     map(seq) do (year, way_names)
#         selected_ways = select_ways(w, way_names)
#         Raster(boolmask(selected_ways; to=dem, shape=:line); name = Symbol("ways_$year"))
#     end |> RasterStack
# end
# plot(way_masks.mus)

# distance_to_roads = map(island_keys, dems, way_masks) do i, dem, ways
#     rast = mask(nearest_distances(ways.primary); with=dem)
#     write(joinpath(distancedir, string(i), "to_primary_roads.tif"), rast)
#     rast = mask(nearest_distances(ways.primary .| ways.secondary); with=dem)
#     write(joinpath(distancedir, string(i), "to_secondary_roads.tif"), rast)
# end


using MakieDraw
using Test
using GLMakie
using GeometryBasics
using GeoInterface
using CSV, DataFrames
using MapRasterization
using GeoJSON
using Colors

includet("raster_common.jl")
includet("map_file_list.jl")

mus_roads_path = joinpath(datadir, "Generated/Roads/mus/mus_roads.tif")
reu_roads_path = joinpath(datadir, "Generated/Roads/reu/reu_roads.tif")

files = define_map_files()

img_path = files.mus.atlas_dutch_period.filename
img_path = files.mus.atlas_18C_land_use.filename
# img_path = files.mus.atlas_19C_land_use.filename
# img_path = files.mus.atlas_1992_land_use.filename
csv_path = splitext(img_path)[1] * ".csv"
points = isfile(csv_path) ? open_warp_points(img_path) : nothing
img = load_image(img_path)
warper = MapRasterization.Warper(; template=dems.mus, point_table=points)
rs = MapRasterization.warp(warper, img; missingval=0)
rs = map(identity, rs)

roadsdir = joinpath(datadir, "Generated/Roads")
roads_dutch_path = joinpath(roadsdir, "roads_dutch.json")
roads_1763_path = joinpath(roadsdir, "roads_1763.json")
roads_added_by_1810_path = joinpath(roadsdir, "roads_added_by_1810.json")
roads_1810_path = joinpath(roadsdir, "roads_1810.json")
roads_1854_path = joinpath(roadsdir, "roads_1854.json")
roads_added_by_1905_path = joinpath(roadsdir, "roads_1905.json")

roads_dutch = GeoJSON.read(read(roads_dutch_path))
roads_1763 = GeoJSON.read(read(roads_1763_path))
roads_added_by_1810 = GeoJSON.read(read(roads_added_by_1810_path))
roads_1810 = GeoInterface.FeatureCollection(
    vcat(collect(GeoInterface.getfeature(roads_1763)), collect(GeoInterface.getfeature(roads_added_by_1810)))
)
roads_1854 = GeoJSON.read(read(roads_1854_path))
roads_added_by_1905 = GeoJSON.read(read(roads_added_by_1905_path))
roads_1905 = GeoInterface.FeatureCollection(
    vcat(collect(GeoInterface.getfeature(roads_1854)), collect(GeoInterface.getfeature(roads_added_by_1905)))
)

const GI = GeoInterface
rod_json = joinpath(datadir, "Roads/rod_ways.geojson")
rod_roads = GI.convert.(Ref(GeometryBasics), filter(GI.geometry.(GeoJSON.read(read(rod_json)))) do x
    GI.trait(x) isa LineStringTrait
end)
Makie.plot(rod_roads)


# GeoJSON.write(roads_1763_path, roads_1763)
# GeoJSON.write(roads_added_by_1810_path, roads_added_by_1810)
# GeoJSON.write(roads_1810_path, roads_1810)
# GeoJSON.write(roads_added_by_1905_path, roads_added_by_1905)
# GeoJSON.write(roads_1854_path, roads_1854)
# GeoJSON.write(roads_dutch_path, line_canvas)

# line_canvas = GeometryCanvas(roads_1763)
# # line_canvas = GeometryCanvas{LineString}()
# Makie.image!(line_canvas.axis, lookup(rs, (X, Y))..., rs)
# GeoJSON.write(line_canvas)

# roads_fcs = (; roads_dutch, roads_1763, roads_1810, roads_1854, roads_1905)
# mus_road_years = [1700, 1763, 1810, 1854, 1905]
# mus_road_series = map(roads_fcs) do fc
#     rebuild(UInt8.(boolmask(fc; to=dems.mus)); missingval=0x00)
# end |> collect |> x -> RasterSeries(x, Ti(mus_road_years))
# plot(mus_road_rasters)
# mus_roads_path = joinpath(datadir, "Generated/Roads/mus/mus_roads.tif")
# write(mus_roads_path, mus_road_rasters)

# # Reunion
# img_path = files.reu.atlas_early_settlement.filename
# img_path = files.reu.atlas_1780_agriculture.filename
# img_path = files.reu.atlas_1815_agriculture.filename
# # img_path = files.reu.atlas_population_1967.filename
# json_roads_path = string(splitext(img_path)[1], "_roads.json")
# aason_roads = GeoJSON.read(read(json_roads_path))
# csv_path = splitext(img_path)[1] * ".csv"
# pts = isfile(csv_path) ? open_warp_points(img_path) : nothing
# img = load_image(img_path)
# rs = MapRasterization.applywarp(img; template=dems.reu, points=pts, missingval=0)
# rs = map(identity, rs)

# # line_canvas = Canvas{LineString}()
# # ax.aspect = AxisAspect(1)
# # draw!(line_canvas, fig, ax)
# line_canvas = GeometryCanvas(json_roads; propertynames=(:built, :closed))
# p = Makie.image!(line_canvas.axis, lookup(rs, (X, Y))..., rs)
# GeoJSON.write(json_roads_path, line_canvas)

# # Generate raster masks
# json_roads = GeoJSON.read(read(json_roads_path))
# reu_road_years = sort(filter(x -> !(x in (" ", "")), union(json_roads.built)))
# df = DataFrame(json_roads)
# reu_road_rasters = map(reu_road_years) do year
#     fc = filter(x -> x.built <= year, df)
#     rebuild(UInt8.(boolmask(fc; to=dems.reu)); missingval=0x00)
# end |> x -> RasterSeries(x, Ti(parse.(Int, reu_road_years)))
# plot(reu_road_rasters, size=(1200, 700))
# reu_roads_path = joinpath(datadir, "Generated/Roads/reu/reu_roads.tif")
# write(reu_roads_path, reu_road_rasters)

# anim = @animate for r in reu_road_rasters
#     plot(dems.reu ./ maximum(skipmissing(dems.reu)))
#     plot!(r)
# end
# gif(anim, "reu_roads.gif", fps=3)
# anim = @animate for r in mus_road_rasters
#     plot(dems.mus ./ maximum(skipmissing(dems.mus)))
#     plot!(r)
# end
# gif(anim, "mus_roads.gif", fps=1)

