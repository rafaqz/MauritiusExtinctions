println("Reading waterways json...")
includet("nearest.jl")

waterways_path = joinpath(datadir, "water.geojson")
waterways_fc = GeoJSON.read(read(waterways_path))
using GeoInterface
println("Converting json to GeometryBasics...")
waterways_rivers = map(waterways_fc) do feature
    geom = GeoInterface.geometry(feature)
    if geom isa GeoJSON.LineString
        GeoInterface.convert(GeometryBasics.LineString, geom)
    else
        missing
    end
end |> skipmissing |> collect
waterways_lakes = map(waterways_fc) do feature
    geom = GeoInterface.geometry(feature)
    if geom isa GeoJSON.Polygon
        GeoInterface.convert(GeometryBasics.Polygon, geom)
    else
        missing
    end
end |> skipmissing |> collect
watermasks = map(dems) do dem 
    rebuild(boolmask(waterways_fc; to=dem); name=:waterways)
end

# distance_to_water = map(island_keys, dems, watermasks) do i, dem, watermask
#     rast = mask(nearest_distances(watermask); with=dem)
#     write(joinpath(distancedir, string(i), "to_water.tif"), rast)
#     rast
# end
distance_to_water = map(island_keys) do k
    Raster(joinpath(distancedir, string(k), "to_water.tif"))
end

# printlnt("Loading national parks json...")
# national_parks = map((mus="mus.geojson", reu="reu.geojson")) do file
#     path = joinpath(datadir, "NationalParks", file) 
#     GeoJSON.read(read(path)).geometry
# end
