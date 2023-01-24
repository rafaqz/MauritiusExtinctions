# using Rasters, Shapefile, DataFrames, Plots, ColorShemes
using Statistics
includet("raster_common.jl")
includet("cost_distance.jl")
includet("map_file_list.jl")
includet("ports.jl")
includet("slope.jl")
includet("roads.jl")

# shp = Shapefile.Table(joinpath(datadir, "Priorisation_actions_de_lutte_-_note_explicative/enjeu_invasion_actions.shp"))
# shp = Shapefile.Table(joinpath(datadir, "Dominique/Past present vegetation shape files/past_vegetation2.shp"))

norder_dir = mkpath(joinpath(outputdir, "Norder"))
norder_stack = mask(replace_missing(RasterStack(norder_dir)[Band(1)]); with=dems.mus)
lc = RasterStack(values(norder_stack)[4:7]...)
plot(lc; size=(1200, 1000))
savefig("landcover.png")
plot(dems.mus ./ maximum(skipmissing(dems.mus)); size=(1200, 1000))
plot!(rebuild(lc[:lc_1773] ./ 2; missingval=0.5); opacity=0.4, title="clearing 1773 with elevation")
savefig("lc_1773_dem.png")
plot!(warped_vegetation; color=:red, linecolor=nothing, title="clearing 1773 with current natives")
savefig("lc_1773_dem_natives.png")

# Generate habitat types from rainfall
# following Strahm

# upland = (norder_stack[:rainfall] .> 2500) .& (dems.mus .> 365)
r = norder_stack[:rainfall] .* u"mm"
e = elevation.mus
palm_savannah = (e .<= 365u"m") .& (r .< 1000u"mm")
upland = (e .> 365u"m") .& (r .> 2500u"mm") # Or is it 3048mm ??
climax_forest = (e .> 365u"m") .& (r .> 3175u"mm") .& (r .< 3556u"mm")
# Mt Cocotte
mossy_forest = (e .> 600"m") .& (r .> 4445"mm")
# phillipa = in.(e, Ref(610u"m" .. 670u"m")) .&  in.(r, Ref(u"m" .. 670u"m"))
sideroxylon = r .> 4400u"mm"
lowland = .!(palm_savannah .| upland)
habitat = RasterStack((; palm_savannah, lowland, upland, sideroxylon, climax))
plot(habitat)
# and Vaughan and Wiehe
palm_savannah=norder_stack[:rainfall] .< 1000
upland = (norder_stack[:rainfall] .> 2500) .& (dems.mus .> 300)
lowland = .!(palm_savannah .| upland)

locations = (; 
    maccabe = (57.443233, -20.392181), 
)
r[map(Contains, locations.maccabe)...]
plot(r)


soiltypespath = joinpath(datadir, "Norder/K factor/SoilK.shp")
soiltypenames = (Symbol.(replace.(Shapefile.Table(soiltypespath).Soil_Group, Ref(" " => "_")))...,)
soilnums = NamedTuple{soiltypenames}((1:14...,))
soilmasks = map(soilnums) do v
    norder_stack[:soiltypes] .== v
end |> RasterStack
plot(soilmasks; c=:viridis)
plot(norder_stack[:soiltypes])
plot(soilmasks; size=(2000, 1000))

# lakesraster = Raster("warpedlakes.tif")[Band(1)]
# elevationraster = Raster("warpedelevation.tif")[Band(1)]
plot(dems.mus ./ maximum(skipmissing(dems.mus)))
plot!(norder_stack[:lc_1935]; c=:viridis, opacity=0.4)
# elevationraster = replace(read(elevationraster), elevationraster[10, 10, 1] => -Inf32)

port_timelines = (
    mus=[1600, 1708],
    reu=[1600, 1886],
)
# distance_stacks = map(island_keys, port_timelines) do i, pti
#     st = read(RasterStack(joinpath(distancedir, string(i)))[Band(1)])
#     to_major_port = Raster(cat(st[:to_major1_ports], st[:to_major2_ports]; dims=Ti(pti)); name=:to_major_ports)
#     return RasterStack(st[(:to_coast, :to_minor_ports, :to_primary_roads,:to_secondary_roads,:to_water)]..., to_major_port) 
# end
# plot(distance_stacks.mus[Ti(2)])
# plot(distance_stacks.mus[:to_secondary_roads])
# plot(distance_stacks.mus)
# plot(distance_stacks.mus[:to_primary_roads])

# Slope
slope_stacks = map(dems) do dem
    slopeaspect(dem, FD3Linear(); cellsize=111.0)
end;
plot(slope_stacks.mus)

files = get_map_files()
slices = make_raster_slices(masks)

function _costs(dems, ports)
    map(dems, ports) do costs, ports
        origins = zeros(Int, dims(costs))
        # map(ports[(:major1, :major2)]) do ps
        map(ports) do ps
            map(ps) do p
                origins[_contains(p)...] = 1
            end
        end
        costfunc = CombinedCost((dem=SlopeCost(slopefactor=-3.5, distfactor=0.05), roads=meancost), *)
        cost_distance(costfunc; origins, costs, cellsize=90)
    end
end

travel_origins = map(dems, ports) do d, p 
    o = map(_ -> Inf * u"hr", d)
    # foreach(dems, travel_origins) do d, o
    # end
    for port in p.major
        o[map(Near, reverse(port))...] = 0.0u"hr"
    end
    for port in p.minor
        o[map(Near, reverse(port))...] = 1.0u"hr"
    end
    o
end

# A history of woods and forests of Mauritius:
# 1 league = 4.83 kilometers per day in forest
# forst_walking_speed = 

resistance = map(
    (mus=slices.mus.timelines.cleared.cleared_1810, reu=slices.reu.timelines.cleared.cleared_1815),
    (mus=way_masks.mus[:ways_1813], reu=way_masks.reu[:ways_1804])) do cleared, ways
    broadcast(cleared, ways) do c, w
        if w
            relative_movement_speed[:dirt_road]
        elseif c
            relative_movement_speed[:cleared_land]
        else
            relative_movement_speed[:forest]
        end
    end
end
elevation = map(dem -> dem .* u"m", dems)

resistance = map(values(slices.mus.timelines.cleared), values(slices.mus.timelines.abandonned)) do cleared, abandonned
    broadcast(cleared, abandonned) do c, a
        if a
            relative_movement_speed_2[:cleared_land]
        elseif c
            relative_movement_speed_2[:cleared_land]
        else
            relative_movement_speed_2[:dense_forest]
        end
    end
end |> NamedTuple{keys(slices.mus.timelines.cleared)}
@time travel_times = map(resistance) do r
    mus_travel_time_map = cost_distance(travel_origins.mus, elevation.mus, r; cellsize=step(lookup(dems.mus, X)) * 111.0u"km")
end
plot(first(travel_times))

anim = @animate for (k, A) in pairs(travel_times)
    plot(A; legend=false, title=string("travel time ", string(k)[end-3:end]), clims=(0, 48))
end
gif(anim, "travel_times.gif", fps=0.8)

@time travel_times = map(travel_origins, elevation.mus, resistance) do o, e, r
    mus_travel_time_map = cost_distance(o, e, r; cellsize=step(lookup(dems.mus, X)) * 111.0u"km")
end
plot(map(t -> plot(t; clims=(0, 40)), travel_times)...)

plot(plot(travel_times.reu), plot(slices.reu.timelines.cleared.cleared_1815), Plots.contourf(dems.reu; levels=[0:500:3000...]))
reu_dist_cleared = replace_missing(travel_times.reu) .* rebuild(replace(slices.reu.timelines.cleared.cleared_1815, 0 => missing); missingval=missing)
mus_dist_cleared = replace_missing(travel_times.mus) .* rebuild(replace(slices.mus.timelines.cleared.cleared_1810, 0 => missing); missingval=missing)
mean(skipmissing(mus_dist_cleared))
mean(skipmissing(reu_dist_cleared))
plot(mus_dist_cleared; clims=(0, 10))
plot(reu_dist_cleared; clims=(0, 10))
count(_ -> true, skipmissing(mus_dist_cleared))
count(_ -> true, skipmissing(reu_dist_cleared))
size(mus_dist_cleared)
size(reu_dist_cleared)
plot(slope.mus)
r

slices.reu.timelines.cleared.cleared_1815 |> sum
slices.mus.timelines.cleared.cleared_1810 |> sum

# Vegetation maps from "Lost Land of the Dodo"
lostland_stacks = map(namedkeys(lostland_image_classes), lostland_image_classes) do i, rasters
    read(RasterStack(joinpath(outputdir, "LostLand", string(i))))
end 
plot(lostland_stacks.reu)
plot(lostland_stacks.mus)

# Vegetation classes
lostland_mask_stacks = map(lostland_stacks, lostland_image_classes) do stack, classes
    map(stack, classes) do r, c
        keys = map(x -> Symbol(replace(x, " " => "_", "-" => "_")), c) |> values
        map((1:length(c)...,)) do id
            classify(r, UInt8(id) => true; others=false, missingval=missing)
        end |> NamedTuple{keys} |> RasterStack
    end |> NamedTuple{keys(stack)}
end

# Cumulative deforestation means deforestation
# by the end of the period. So we rename.
# deforestation_phases = (
#     mus = (
#        by_1807=:deforested_before_1807, 
#        by_1835=:deforested_1807_1835,
#        by_1910=:deforested_1835_1910,
#        by_1947=:deforested_1910_1947,
#        by_1970=:deforested_1947_1970,
#        by_2010=:deforested_since_1970
#     ),
#     reu = (
#        by_1700=:C17, 
#        by_1800=:C18,
#        by_1900=:C19,
#        by_2000=:C20,
#     ),
# )
# deforestation = map(lostland_mask_stacks, deforestation_phases) do stack, phases
#     phase_stack = stack.phase[values(phases)]
#     reduce(NamedTuple(phase_stack); init=()) do acc, A
#         acc === () ? (A,) : (acc..., last(acc) .| A) 
#     end |> xs -> RasterStack((missingmask(first(xs)), xs...); keys=(:by_1600, keys(phases)...))
# end

# Landcover
lc_dir = joinpath(datadir, "Landcover/")
lc_names = ( :No_Data,
  :Continuous_urban,
  :Discontinuous_urban,
  :Forest,
  :Shrub_vegetation,
  :Herbaceaous_vegetation,
  :Mangrove,
  :Barren_land,
  :Water,
  :Sugarcane,
  :Pasture,
  :UnusedIndex,
  :Other_cropland,
)
lc_categories = NamedTuple{lc_names}((Int32.(0:12)...,))

# Read from tif
lc_rasterized = map(island_keys) do island
    path = joinpath(lc_dir, "$(island)_landcover.tif")
    rast = Raster(path; name=:modern_landcover)[Band(1)]
    # Get rid of Water and NoData
    replace_missing(replace_missing(rast, lc_categories.Water), lc_categories.No_Data)
end
plot(lc_rasterized.mus)

# Masks for each land cover
lc_masks = map(lc_rasterized) do rast
    masks = map(lc_categories) do v
        mask(rast .== v; with=rast, missingval=missing)
    end
    RasterStack(masks)
end
lc_2017 = map(lc_masks) do m 
    stack = RasterStack((;
        urban=m[:Continuous_urban] .| m[:Discontinuous_urban],
        forest=m[:Shrub_vegetation] .| m[:Forest],
        cleared=m[:Pasture] .| m[:Sugarcane] .| m[:Other_cropland],
    ))
    other = .!(reduce((a, b) -> a .| b, values(stack)))
    RasterStack((; stack..., other))
end 
plot(lc_2017.mus)

plot(lc_2017.mus; size=(2000, 1000))
plot(lc_2017.mus[:forest]; size=(2000, 1000))
p = plot(lc_2017.mus[:forest])
plot!(soilmasks[:Lithosol] .| soilmasks[:Lithosol]; opacity=0.5)
p = plot(soilmasks[:Lithosol] .| soilmasks[:Lithosol]; opacity=0.8)
for i in 1:3
    class = filter(row -> row.category == i, collect(warped_vegetation))
    c = first(class).color
    color = RGB(c[:r], c[:g], c[:b])
    plot!(class; alpha=0, color, fillalpha=0.5)
end
display(p)
# plot(lc_masks.mus; c=:viridis)





using Colors, GeometryBasics
using GLMakie
using CairoMakie
using XML
xml = XML.Document("/home/raf/PhD/Mascarenes/Data/Selected/Mauritius/forest.svg");

allpaths = filter(children(xml.root)) do e
    e isa Element && tag(e) == "path"
end

color_strings = map(allpaths) do x
    hasproperty(x, :fill) ? x.fill : "none"
end |> union

category_paths = map(color_strings) do color_string
    color = if color_string == "none"
        RGB(0.0, 0.0, 0.0)
    else
        RGB(map(s -> parse(Float64, s) / 100, split(color_string[5:end-2], "%, "))...)
    end
    paths = filter(x -> x.fill == color_string, allpaths)
    path_strings = map(p -> p.d, paths)
    linestrings = Vector{Tuple{Float64,Float64}}[]
    geoms = map(path_strings) do p
        points = map(filter(!isempty, split(p, "L "))) do c
            points = map(x -> parse(Float64, x), filter(!in(("", "Z", "M")), split(c, " ")))
            Point2(points[1], points[2])
        end |> x -> filter(!isempty, x)
        if occursin("Z", p)
            Polygon(points)
        else
            LineString(points)
        end
    end |> x -> filter(!isempty, x)
    (; color, geoms)
end

colors = getproperty.(category_paths, :color)
fig = GLMakie.Figure(resolution = (1800,1600))
fig[1,1] = ax = Axis(fig)
ax.xreversed = false
ax.yreversed = false
c = category_paths[2]
CairoMakie.lines!(ax, c.geoms; color=c.color)
for c in category_paths[4:7]
    if typeof(first(c.geoms)) <: Polygon
        GLMakie.poly!(ax, c.geoms; color=c.color)
    else
        GLMakie.lines!(ax, c.geoms; color=c.color)
    end
end
display(fig)

save("plot.png", fig, px_per_unit = 4)



methodswith(typeof(first(eachelement(r))))
elements(r)
namespace(r)
