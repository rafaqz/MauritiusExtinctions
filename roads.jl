using Dates, DataStructures
using DimensionalData.LookupArrays
using GeometryBasics
using GeoInterface
using GeoInterfaceRecipes
GeoInterfaceRecipes.@enable_geo_plots GeoJSON.Feature

includet("raster_common.jl")
includet("waynames.jl")

# To roads
function namelike(names, feature)
    haskey(feature.properties, "name") || return false
    return any(occursin(feature.properties["name"]), names)
end
function roadin(feature, roadlist)
    hasproperty(feature, :highway) || return false
    feature.highway in ("trunk", "primary", "secondary", "tertiary", "unclassified", "track", "residential") || return false
    return hasproperty(feature, :ref) && getproperty(feature, :ref) in roadlist ||
        hasproperty(feature, :name) && getproperty(feature, :name) in roadlist ||
        hasproperty(feature, Symbol("@id")) && getproperty(feature, Symbol("@id")) in roadlist
end
select_ways(ways, way_names) =
    GeoJSON.FeatureCollection([f for f in ways if roadin(f, way_names)])
function plot_ways(ways, way_names, dem, title="")
    selected_ways = select_ways(ways, way_names)
    plot(dem ./ maximum(skipmissing(dem)), legend=false)
    # plot((dems.mus ./ maximum(skipmissing(dems.mus)))[X=57.4..58.6, Y=(-20.35)..(-20.15)])
    # plot!(watermasks.mus; color=:blue, legend=false) 
    # plot!(deforestation.mus[:by_1835]; c=:viridis, opacity=0.4) 
    # plot!(selected_ways; color=:red, title, legend=false)
    plot!(selected_ways; title, legend=false)
end
function rasterize_ways(ways, way_names; to)
    selected_ways = select_ways(ways, way_names)
    rasterize(selected_ways; to, fill=true) 
end

ways = map(island_keys) do ik
    json_path = joinpath(datadir, "Roads", "$(ik)_ways.geojson")
    GeoJSON.read(read(json_path))
end

way_names_sequences = (; 
    mus = (
        1667=>mus_way_names_1667,
        1725=>mus_way_names_1725,
        1764=>mus_way_names_1764,
        1813=>mus_way_names_1813,
        1880=>mus_way_names_1880,
    ), 
    reu = (
        1752=>reu_way_names_1752,
        1804=>reu_way_names_1804,
    )
)

selected_ways = select_ways(ways.mus, way_names_sequences.mus[4][2])

road_lines = map(ways) do island_ways
    map(island_ways) do way
        geom = GeoInterface.geometry(way)
        if geom isa GeoJSON.LineString
            GeoInterface.convert(GeometryBasics.LineString, geom)
        else
            missing
        end
    end |> skipmissing |> collect
end 

way_masks = map(ways, way_names_sequences, dems) do w, seq, dem
    map(seq) do (year, way_names)
        selected_ways = select_ways(w, way_names)
        Raster(boolmask(selected_ways; to=dem, shape=:line); name = Symbol("ways_$year"))
    end |> RasterStack
end
plot(way_masks.mus)

# distance_to_roads = map(island_keys, dems, way_masks) do i, dem, ways
#     rast = mask(nearest_distances(ways.primary); with=dem)
#     write(joinpath(distancedir, string(i), "to_primary_roads.tif"), rast)
#     rast = mask(nearest_distances(ways.primary .| ways.secondary); with=dem)
#     write(joinpath(distancedir, string(i), "to_secondary_roads.tif"), rast)
# end
