# using Rasters, Shapefile, DataFrames, Plots, ColorShemes
using Statistics, Chain
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
# Plots.plot(norder_stack[:lc_1935])
# Plots.plot!(warped_vegetation)
# plot(lc; size=(1200, 1000))
# savefig("landcover.png")
# plot(dems.mus ./ maximum(skipmissing(dems.mus)); size=(1200, 1000))
# plot!(rebuild(lc[:lc_1773] ./ 2; missingval=0.5); opacity=0.4, title="clearing 1773 with elevation")
# savefig("lc_1773_dem.png")
# plot!(warped_vegetation; color=:red, linecolor=nothing, title="clearing 1773 with current natives")
# savefig("lc_1773_dem_natives.png")

# Generate habitat types from rainfall
# following Strahm

# upland = (norder_stack[:rainfall] .> 2500) .& (dems.mus .> 365)
r = norder_stack[:rainfall] .* u"mm"
e = elevation.mus
palm_savannah = (e .<= 365u"m") .& (r .< 1000u"mm")
upland = (e .> 365u"m") .& (r .> 2500u"mm") # Or is it 3048mm ??
climax_forest = (e .> 365u"m") .& (r .> 3175u"mm") .& (r .< 3556u"mm")
# Mt Cocotte
mossy_forest = (e .> 600u"m") .& (r .> 4445u"mm")
# phillipa = in.(e, Ref(610u"m" .. 670u"m")) .&  in.(r, Ref(u"m" .. 670u"m"))
sideroxylon = r .> 4400u"mm"
lowland = .!(palm_savannah .| upland)
habitat = RasterStack((; palm_savannah, lowland, upland, sideroxylon, climax_forest))
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
plot(soilmasks)
lith = skipmissing(soilmasks.Lithosol)
reg = skipmissing(soilmasks.Regosol)
map(soilmasks) do A
    sum(skipmissing(A)) / count(x -> true, skipmissing(A))
end

soilslopestack = merge(soilmasks, slope_stacks.mus, (; elevation=dems.mus))
soildf = DataFrame(soilslopestack)
soil_fit = map(slices.mus.timelines.cleared) do a
    soildf[!, :x] = vec(parent(a))
    forula = FormulaTerm(Term(:x), map(Term, keys(soilslopestack)))
    model = lm(forula, soildf)
end;
soil_fit[12] |> r2
plot(soilmasks.Dark_Magnesium_Clay)
plot(RasterStack(slices.mus.timelines.cleared[3]))
plot(soilmasks.Lithosol)
plot(mask(slices.mus.timelines.cleared[3]; with=masks.mus))
s = map(1:13) do i
    soil_cleared = map(soilmasks) do A
        m = mask(A; with=rebuild(slices.mus.timelines.cleared[i]; missingval=false))
        sum(skipmissing(m)) / sum(skipmissing(A))
    end |> pairs
end
s[6]
plot(soil_cleared)
# plot(soilmasks; c=:viridis)
# plot(norder_stack[:soiltypes])
# plot(soilmasks; size=(2000, 1000))

# lakesraster = Raster("warpedlakes.tif")[Band(1)]
# elevationraster = Raster("warpedelevation.tif")[Band(1)]
# plot(dems.mus ./ maximum(skipmissing(dems.mus)))
# plot!(norder_stack[:lc_1935]; c=:viridis, opacity=0.4)
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
lc_names = (
  :No_Data,
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
lc_2017_category_groups = (
    forest_or_abandoned = (:Forest, :Mangrove, :Shrub_vegetation, :Barren_land),
    urban = (:Continuous_urban, :Discontinuous_urban),
    cleared = (:Sugarcane, :Other_cropland, :Pasture),
    uncertain = (:Herbaceaous_vegetation,),
)

lc_2017_categories = NamedTuple{lc_names}((Int32.(0:12)...,))

# Read from tif
lc_2017_rasterized = map(island_keys) do island
    path = joinpath(lc_dir, "$(island)_landcover.tif")
    rast = Raster(path; name=:modern_landcover)[Band(1)]
    # Get rid of Water and NoData
    replace_missing(replace_missing(rast, lc_2017_categories.Water), lc_2017_categories.No_Data)
end
plot(lc_2017_rasterized.mus)

# Masks for each land cover
lc_2017_masks = map(lc_2017_rasterized) do rast
    masks = map(lc_2017_categories) do v
        mask(rast .== v; with=rast, missingval=false)
    end
    RasterStack(masks)
end
lc_2017 = map(lc_2017_masks) do m 
    stack = map(lc_2017_category_groups) do group
        slices = m[group]
        reduce((acc, x) -> acc .| x, values(slices))
    end |> RasterStack
end 
plot(lc_2017.reu)


land_use_2002_shpfile = joinpath(datadir, "Claudia/Demo/GIS WILD LIFE FOUNDATION/SHAPE FILE/land_use_WGS_region.shp")
lu_2002 = Shapefile.Table(land_use_2002_shpfile)
lu_2002df = DataFrame(lu_2002)
lu_2002_categories = map(enumerate(union(skipmissing(lu_2002df.LAND_USE)))) do (i, cat)
    Symbol(replace(cat, " " => "_")) => (i => cat)
end |> NamedTuple
lu_2002_categories |> pairs
dfs = map(lu_2002_categories) do cat
    df = filter(lu_2002df) do row
        !ismissing(row.LAND_USE) && row.LAND_USE == cat[2]
    end
    df[!, :category] .= cat[1]
    df
end
lu_2002df_categorised = vcat(dfs...)
lu_2002_rast = rasterize(lu_2002df_categorised; to=dems.mus, fill=:category)
lu_2002_swamp = lu_2002_rast .== lu_2002_categories.Marsh_or_Swamp[1]
plot(lu_2002_swamp)
lu_2002_agriculture = lu_2002_rast .!= lu_2002_categories.Marsh_or_Swamp[1]
plot(lu_2002_agriculture)

# 2002 agriculture
lu_2002_agriculture |> plot
# 2002 Agriculture cleared in 1992
slices.mus.timelines.cleared.cleared_1992 .& lu_2002_agriculture |> plot
slices.mus.timelines.cleared.cleared_1992 .& .!(lu_2002_agriculture) |> plot
.!(slices.mus.timelines.cleared.cleared_1992) .& lu_2002_agriculture |> plot
# Cleared areas not in agriculture in 2002
(slices.mus.timelines.lc.lc_1992 .== lc_categories.cleared) .& .!(boolmask(lu_2002_agriculture .| lu_2002_swamp)) |> plot

includet("svgs.jl")
mus_native_2017_lc = mask(lc_2017_rasterized.mus; with=mus_native_veg_rast)
countcats(mus_native_2017_lc, lc_2017_categories) |> pairs
mus_invasives_2017_mask = (.!(mus_native_veg_mask) .& (lc_2017.mus.forest_or_abandoned .== 1))
plot(mus_invasives_2017_mask)
mus_forestry_1992 = rebuild(mask(slices.mus.timelines.forestry.forestry_1992 .* .!(mus_native_veg_mask); with=masks.mus); name=:forestry)
plot(mus_forestry_1992)
mus_uncertain_inasives_2017_mask = (.!(mus_native_veg_mask) .& (lc_2017.mus.uncertain.== 1))
mus_all_invasives_2017 = (mus_invasives_2017_mask .+ mus_uncertain_inasives_2017_mask .* 0.5) .* .!(forestry_1992)
mus_invasive_density = @chain begin
    map(mus_all_invasives_2017, mus_native_density) do i, n
        !ismissing(n) && n > 0 ? (1 - n) : i
    end
    rebuild(_; name=:invasive_density) 
    mask(_; with=dems.mus, missingval=missing)
end
plot(RasterStack(mus_native_density, mus_invasive_density, mus_forestry_1992);
     size=(1000, 1000), clims=(0, 1),
)
savefig("mauritius_vegtetation.png")


plot(lc_2017.mus; size=(2000, 1000))
plot(lc_2017.mus.forest]; size=(2000, 1000))
p = plot(lc_2017.mus.forest_or_abandoned)
plot!(soilmasks.Lithosol; opacity=0.5)
p = plot(soilmasks.Lithosol .| soilmasks.Regosol; opacity=0.8)
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
