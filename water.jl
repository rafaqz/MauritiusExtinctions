println("Reading waterways json...")
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

# printlnt("Loading national parks json...")
# national_parks = map((mus="mus.geojson", reu="reu.geojson")) do file
#     path = joinpath(datadir, "NationalParks", file) 
#     GeoJSON.read(read(path)).geometry
# end
