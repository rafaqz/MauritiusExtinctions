
batcaves = CSV.read("/home/raf/Downloads/cave_locations_Ashmi.csv", DataFrame)

points = collect(zip(batcaves.long, .-(batcaves.lat)))
data = extract(NamedTuple.(complete.timeline[Ti(1600..2017)]), points)
dfs = map(data, batcaves.cave) do ((long, lat), vals), name
    map(vals, lookup(vals, Ti)) do val, time
        merge((; name, time, long, lat), NamedTuple(val))
    end |> DataFrame
end
locs = vcat(dfs...)
p = Rasters.rplot(striped_complete[Ti=(At(1965))])
Makie.scatter!(p.axis, points)
save("landcover_with_caves.png", p)
CSV.write("cave_landcover_history.csv", locs)

