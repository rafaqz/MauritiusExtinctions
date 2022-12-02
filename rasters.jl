# using Rasters, Shapefile, DataFrames, Plots, ColorShemes
includet("raster_common.jl")
includet("cost_distance.jl")

# shp = Shapefile.Table(joinpath(datadir, "Priorisation_actions_de_lutte_-_note_explicative/enjeu_invasion_actions.shp"))
# shp = Shapefile.Table(joinpath(datadir, "Dominique/Past present vegetation shape files/past_vegetation2.shp"))
gc = df[!, :geometry]
colors = getindex.(Ref(ColorSchemes.viridis), gc)
p = plot()
length()
for (i, shp) in enumerate(df[!, :geometry])
    plot!(p, df[!, :geometry][1]; col=colors[i][1])
end
display(p)

minimum(gc)
plot(shp; )

df = DataFrame(shp)
names(df)
union(df[!, :Invasion])

norder_dir = mkpath(joinpath(outputdir, "Norder"))
norder_stack = mask(replace_missing(RasterStack(norder_dir)[Band(1)]); with=dems.mus)

soiltypespath = joinpath(datadir, "Norder/K factor/SoilK.shp")
soiltypenames = (Symbol.(replace.(Shapefile.Table(soiltypespath).Soil_Group, Ref(" " => "_")))...,)
soilnums = NamedTuple{soiltypenames}((1:14...,))
soilmasks = map(soilnums) do v
    norder_stack[:soiltypes] .== v
end |> RasterStack
plot(soilmasks; c=:viridis)
plot(norder_stack[:soiltypes])

# lakesraster = Raster("warpedlakes.tif")[Band(1)]
# elevationraster = Raster("warpedelevation.tif")[Band(1)]
plot(dems.mus ./ maximum(skipmissing(dems.mus)))
plot!(norder_stack[:lc_1835]; c=:viridis, opacity=0.4)
# elevationraster = replace(read(elevationraster), elevationraster[10, 10, 1] => -Inf32)

port_timelines = (
    mus=[1600, 1708],
    reu=[1600, 1886],
)
distance_stacks = map(island_keys, port_timelines) do i, pti
    st = read(RasterStack(joinpath(distancedir, string(i)))[Band(1)])
    to_major_port = Raster(cat(st[:to_major1_ports], st[:to_major2_ports]; dims=Ti(pti)); name=:to_major_ports)
    return RasterStack(st[(:to_coast, :to_minor_ports, :to_primary_roads,:to_secondary_roads,:to_water)]..., to_major_port) 
end
distance_stacks.mus[:to_major_ports]
plot(distance_stacks.mus[:to_secondary_roads])


# Slope
slope_stacks = map(dems) do dem
    slopeaspect(dem, FD3Linear())
end
plot(slope_stacks.reu)

# Gives a weird answer
# slopecost(slope) = 6 / exp(-3.5 * (slope + 0.05))
# slopecost(0.1)
using BenchmarkTools
using ProfileView

_contains((y, x)) = Y(Contains(y)), X(Contains(x))
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

includet("ports.jl")
ag = 2
@time costs = _costs(dems, ports)
@time agcosts = _costs(agdems, ports)
agcosts2 = map(d -> Rasters.aggregate(mean, d, ag), costs)
costmean = maximum((maximum(skipmissing(agcosts.reu)), maximum(skipmissing(agcosts2.reu))))
plot((agcosts.reu .- agcosts2.reu) / costmean; )

c_d.mus
clims = (0, 3000)
color = :magma
plot(plot(c_d.reu; clims, color, title="mus cost"), plot(c_d.mus; clims, color, title="reu cost"); size=(1000,1000))

savefig("costdistance_allports.png")

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

typeof(lostland_image_classes)

# Cumulative deforestation means deforestation
# by the end of the period. So we rename.
deforestation_phases = (
    mus = (
       by_1807=:deforested_before_1807, 
       by_1835=:deforested_1807_1835,
       by_1910=:deforested_1835_1910,
       by_1947=:deforested_1910_1947,
       by_1970=:deforested_1947_1970,
       by_2010=:deforested_since_1970
    ),
    reu = (
       by_1700=:C17, 
       by_1800=:C18,
       by_1900=:C19,
       by_2000=:C20,
    ),
)
deforestation = map(lostland_mask_stacks, deforestation_phases) do stack, phases
    phase_stack = stack.phase[values(phases)]
    reduce(NamedTuple(phase_stack); init=()) do acc, A
        acc === () ? (A,) : (acc..., last(acc) .| A) 
    end |> xs -> RasterStack((missingmask(first(xs)), xs...); keys=(:by_1600, keys(phases)...))
end

plot(deforestation.mus; c=:viridis)
plot(deforestation.reu; c=:viridis)
plot(distance_stacks.mus)

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

# plot(lc_masks.mus; c=:viridis)

nothing



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
