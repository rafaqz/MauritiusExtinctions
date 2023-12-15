using GeoInterfaceRecipes
using Shapefile, DataFrames, GBIF2, Plots
using GeoJSON
using CSV
using MapRasterization
using Colors
using FileIO
using ImageIO
using GLMakie
using GeometryBasics
using JSON3
using ColorSchemes
using GeoInterfaceRecipes
using Statistics
GeoInterfaceRecipes.@enable_geo_plots Polygon

includet("raster_common.jl")
# includet("roads.jl")
includet("water.jl")
includet("map_file_list.jl")
# includet("svgs.jl")

ch_original_veg_path= "/home/raf/PhD/Mascarenes/Data/Selected/Mauritius/Undigitised/page33_mauritius_vegetation_colored.tif"
ch_original_veg_raster = Raster(ch_original_veg_path; missingval=-9223372036854775808)
Makie.heatmap(ch_original_veg_raster)
ch_original_veg_json = "/home/raf/PhD/Mascarenes/Data/Selected/Mauritius/Undigitised/page33_mauritius_vegetation_colored.json"
layer_names = JSON3.read(ch_original_veg_json).settings.category_name


#############################################################################3
# Mauritius vegetation

mus_native_veg_df = DataFrame(mus_native_veg_poly)
mus_native_veg_df.geometry = GeoInterface.convert.(GeometryBasics.Polygon, mus_native_veg_df.geometry)
mus_native_veg_rast = rasterize(maximum, mus_native_veg_df; to=dems.mus, fill=:category, boundary=:touches)
mus_native_veg_tif_path = "/home/raf/PhD/Mascarenes/Data/Generated/mus_native_veg.tif"
Rasters.rplot(mus_native_veg_rast)
write(mus_native_veg_tif_path, replace_missing(mus_native_veg_rast, 0); force=true)
mus_native_density = rasterize(maximum, mus_native_veg_df; to=dems.mus, fill=:category, boundary=:touches)
mus_native_veg_mask = boolmask(mus_native_veg_rast)
grade_fractions = (0.2, 0.5, 0.7)
mus_native_density = mask(rebuild(map(mus_native_veg_rast) do x
    ismissing(x) || x == 0 ? 0.0 : grade_fractions[x]
end; name=:native_density_1999); with=dems.mus, missingval=missing)
Plots.plot(mus_native_density)

set_theme!(palette = (color = ColorSchemes.tab20, patchcolor=RGBA.(ColorSchemes.tab20, 0.5)))
load_image(img_path::String) = RGB{Float64}.(load(img_path) |> rotr90)
mus_vegmap_path = "/home/raf/PhD/Mascarenes/Data/Selected/Mauritius/Undigitised/VaugnWiehe_Vegetation.png"
img = load_image(mus_vegmap_path)
# output = MapRasterization.select_category_shapes(Polygon, img)
# output = MapRasterization.select_category_shapes(Polaygon, output, img)
output_path = splitext(mus_vegmap_path)[1] * ".json"
output = JSON3.read(read(output_path), MapRasterization.CategoryShapes{Polygon})
output = MapRasterization.select_category_shapes(output, img)
vw_veg_json = "/home/raf/PhD/Mascarenes/Data/Selected/Mauritius/Undigitised/VaugnWiehe_Vegetation.geojson"
vw_fc = map(eachrow(DataFrame(output))) do row
    GeoInterface.Feature(row.geometry; properties=(category=row.category, name=row.name))
end |> GeoInterface.FeatureCollection
GeoJSON.write(vw_veg_json, vw_fc)
vw_fc = GeoJSON.read(vw_veg_json)
vw_df = DataFrame(vw_fc)

guide = (borders.mus, waterways_rivers, waterways_lakes, mus_native_veg_poly...);
vw_ecotypes = warp_to_raster(mus_vegmap_path, dems.mus .* 0; 
    object_type=MapRasterization.CategoryShapes{Polygon}, edit=false, guide
) |> DataFrame

vw_ecotypes_categories = Base.union(vw_ecotypes.name)

ecotype_group_names = (
    native = (
        "Lowland indigenous thicket",
        "Dry forest",
        "Transitional forest",
        "Upland forest",
        "Marshes",
        "Mosay forest",
        "Sideroxylon scrub",
        "Pole forest",
    ),
    lost = (
        "Marshes",
        "Mosay forest",
        "Sideroxylon scrub",
    ),
    invasive = (
        "Cordia scrub",
        "Savannah thorn scrub",
        "Leucaena thicket",
        "Albizzia forest",
        "Psidium scrub",
        "Ravenala thicket",
        "Upland thicket",
        "Grassland",
        "Litsea thicket",
    ),
    plantation = (
        "Plantation",
        "Casuarina plantation",
    ),
    psidium = (
        "Psidium scrub",
    ),
)

ecotype_names = Base.union(vw_ecotypes.name)
map(ecotype_group_names) do names
    map(names) do name
        name in Base.union(vw_ecotypes.name)
    end |> all
end

ecotype_groups = map(ecotype_group_names) do names
    filter(r -> r.name in names, vw_ecotypes)
end
ecotype_group_masks = map(ecotype_groups) do df
    mask(boolmask(df; to=dems.mus); with=dems.mus)
end |> RasterStack
Plots.plot(ecotype_group_masks)
ecotype_rasters = map(ecotype_groups) do df
    mask(rasterize(last, df; to=dems.mus, fill=:category); with=dems.mus)
end |> RasterStack
Plots.plot(ecotype_rasters; 
     color=palette(:batlow,  UnitRange(extrema(skipmissing(ecotype_rasters.invasive))...)),
)
Plots.plot(ecotype_rasters.invasive; color=palette(:thermal, UnitRange(extrema(skipmissing(ecotype_rasters.invasive))...)))
# masked_ecotype_rasters = map(ecotype_rasters) do A
# end
m = boolmask(slices.mus.timelines.cleared.cleared_1905; missingval=0)
Plots.plot(m)

native_1937_alread_cleared = rebuild(replace_missing(ecotype_rasters.native .* m, 0); name="actually_cleared")
Plots.plot(Plots.plot(native_1937_alread_cleared), Plots.plot(ecotype_rasters.native))

# Generate 1937 native density from Vaugn and Whie
mus_native_density_1937 = let m = ecotype_group_masks
    map(m.native, m.invasive, m.plantation) do n, i, p 
        if p 
            0.0
        elseif n & i
            0.5
        elseif i
            0.0
        elseif n
            1.0
        else
            0.0
        end
    end
end |> x -> rebuild(mask(x; with=dems.mus, x, missingval=missing); name=:native_density_1937)
# Correct for poor outlines in the 1937 map by taking the maximum quality of 1937 and 1999. 
# This does *not* correct for where the 1937 area is too large.
# There are few changes, mostly some edges around ferney and tamarin gorge
mus_native_density_1937 .= max.(mus_native_density_1937, mus_native_density)
# 1600 native density is just a mask of 1s
mus_native_density_1600 = rebuild(mus_native_density_1937 .> -1; name=:native_density_1600)
native_density = RasterStack(mus_native_density_1600, mus_native_density_1937, mus_native_density)
Plots.plot(native_density; c=:bamako, layout=(1, 3), size=(1000, 300), clims=(0, 1))
savefig("mus_native_density.png")

#
# GBIF plant locations
#

using CSV, DataFrames, GBIF2
country_codes_2 = (mus=:MU, reu=:RE)

gbif_mauritius_natives = "/home/raf/PhD/Mascarenes/Tables/GBIF/Mauritius_native_plants/taxon.txt"
mus_gbif_native_list = CSV.read(gbif_mauritius_natives, DataFrame)
mus_gbif_native_list.specificEpithet
names(mus_gbif_native_list)
native_scientific_names = map(mus_gbif_native_list.scientificName) do s
    length(split(s)) > 1 ? join(split(s, ' ')[1:2], ' ') : s
end
filter(native_scientific_names) do sf
    occursin("Premna", sf)
end


gbif_plants_csv_path = (
    mus = "/home/raf/PhD/Mascarenes/Tables/mus_gbif_plants.csv",
    reu = "/home/raf/PhD/Mascarenes/Tables/reu_gbif_plants.csv",
)
# mus_gbif_plants = occurrence_search("Tracheophyta"; country=:MU, limit=100000) |> DataFrame
# CSV.write(gbif_plants_csv.mus, mus_gbif_plants)
mus_gbif_plants = CSV.read(gbif_plants_csv_path.mus, DataFrame)
mus_gbif_plants.scientificName
names(mus_gbif_plants)
mus_native_locations = filter(mus_gbif_plants) do p
    p.scientificName in mus_gbif_native_list.scientificName &&
    # p.order in ("Myrtales", "Gentianales", "Ericales", "Rosales")
    !(p.order in ("Poales", "Gentianales", "Malvales")) && # Remove herbs
    !(p.genus in ("Portulaca", "Lactuca", "Striga")) &&
    !(p.species in ("Pandanus utilis", "Premna serratifolia")) && # Naturalised?
    !(p.species in ("Ipomoea pes-caprae", "Scaevola sericea")) && # pantropical dune plants
    !ismissing(p.decimalLatitude)
end
map(union(mus_native_locations.genus)) do genus
    count(==(genus), mus_native_locations.genus) => genus
end |> sort
map(union(mus_native_locations.order)) do order
    count(==(order), mus_native_locations.order) => order
end |> sort

# all_tracheophyta = map(country_codes_2) do country
#     println(country)
#     occurrence_search("Tracheophyta"; country, limit=100000)
# end

using Tyler, MapTiles, TileProviders
# df = mus_gbif_plants
df = mus_native_locations
tuple_points = [p for p in zip(df.decimalLongitude, df.decimalLatitude) if !ismissing(p[1])]
points = map(tuple_points) do pt
    Point2f(MapTiles.project(pt, MapTiles.wgs84, MapTiles.web_mercator))
end
tyler = Tyler.Map(Rect2f(57.0, -21.5, 0.24, 0.15); provider=Google())
# xa = filter(x -> !ismissing(x.genus) && x.genus == "Labourdonnaisia", all_trach_mus)
# p = plot(native_density.native_density_1999; c=:bamako)
# p = Makie.heatmap(lookup(native_density)..., parent(native_density.native_density_1999); transparent=true)
# p = Makie.heatmap(lookup(native_density)..., parent(lc_2017.mus.forest_or_abandoned); transparent=true, alpha=0.2)
Makie.scatter!(tyler.axis, points)#; markersize=1)
Makie.text!(tyler.axis, points; text=df.species, color=:white) 
p.axis.aspect = AxisAspect(1)

# gbif_layer = Leaflet.Layer(GeoInterface.convert.(Val(:GeoJSON), tuple_points); color=:red)
# quality_layers = map((x, c) -> Leaflet.Layer(GeoInterface.convert.(Val(:GeoJSON), x); color=c), veg_quality_1999, (:blue, :green, :yellow))
# ecotypes_layers = map((l, color) -> Leaflet.Layer(GeoInterface.convert.(Val{:GeoJSON}(), l.geometry); color), ecotype_groups[(:native, :invasive)], (native=:orange, invasive=:red))
# layers = [quality_layers..., ecotypes_layers.native, ecotypes_layers.invasive, gbif_layer]
# m = Leaflet.Map(; layers, provider, zoom=3, height=1000, center=[-20.0, 57.0]);
# w = Blink.Window(; body=m)

df = DataFrame(all_plants)
names(DataFrame(all_plants))

CSV.write("/home/raf/PhD/Mascarenes/Tables/reu_gbif_plants.csv", all_plants)
mus_gbif_plants = CSV.read("/home/raf/PhD/Mascarenes/Tables/mus_gbif_plants.csv", DataFrame)
points = collect(zip(mus_gbif_plants.decimalLongitude, mus_gbif_plants.decimalLatitude))
plot(dems.mus)
scatter!(points)
annotate!(points)


#############################################################################3
# Regunion vegetation

reu_past_veg_raster = Raster(joinpath(datadir, "Dominique/Vegetation_Rasters/pastveg3.tif")) |>
    x -> resample(x; to=dems.reu)
reu_present_veg_raster = Raster(joinpath(datadir, "Dominique/Vegetation_Rasters/present_veg4.tif")) |>
    x -> resample(x; to=dems.reu)
plot(reu_past_veg_raster)
plot(reu_present_veg_raster)

# reu_past_veg_path = "/home/raf/PhD/Mascarenes/Data/Dominique/Past present vegetation shape files/past_vegetation2.shp"
# reu_present_veg_path = "/home/raf/PhD/Mascarenes/Data/Dominique/Past present vegetation shape files/present_veg2.shp"
# # Projected in Gauss_Schreiber Transverse Mercator
# # https://proj.org/operations/projections/gstmerc.html
# using GeoFormatTypes, ArchGDAL, GeoDataFrames, OrderedCollections
# sourcecrs = GeoFormatTypes.EPSG(3727) # Reunion 1947
# targetcrs = GeoFormatTypes.EPSG(4326)
# reu_past_veg = GeoDataFrames.read(reu_past_veg_path)
# sort(reu_past_veg, :GRIDCODE)
# reu_present_veg = GeoDataFrames.read(reu_present_veg_path)
# ArchGDAL.reproject(reu_past_veg.geometry, sourcecrs, targetcrs; order=:trad)
# ArchGDAL.reproject(reu_present_veg.geometry, sourcecrs, targetcrs; order=:trad)
# reu_past_veg_raster = rasterize(reu_past_veg; to=dems.reu, fill=:GRIDCODE, missingval=0, boundary=:touches)
# plot(reu_past_veg_raster)
# plot(reu_past_veg.geometry; opacity=0.5)
# plot(reu_present_veg.geometry; opacity=0.5)
# reu_veg_categories = OrderedDict(sort(Base.union(reu_past_veg.NOM) .=> Base.union(reu_past_veg.GRIDCODE); by=last))
# present_GRIDCODE = map(reu_present_veg.NOM) do nom
#     if ismissing(nom) || nom == "no data" || nom == "no_data"
#         0
#     elseif haskey(reu_veg_categories, nom)
#         reu_veg_categories[nom]
#     else
#         reu_veg_categories[nom] = maximum(values(reu_veg_categories)) + 1 
#     end
# end
# reu_veg_categories
# reu_present_veg.GRIDCODE = present_GRIDCODE
# reu_present_veg_raster = rasterize(reu_present_veg; to=dems.reu, fill=:GRIDCODE, missingval=0, boundary=:touches)
# plot(plot(reu_past_veg_raster; title="past veg", clims=(0, 21)) ,
#      plot(reu_present_veg_raster; title="present veg", clims=(0, 21)); size=(1000, 400))
# savefig("reunion_vegetation.png")

# p = plot(reu_past_veg_raster .* 0; title="past veg")
# for i in 1:maximum(values(reu_veg_categories))
#     layer = filter(x -> x.GRIDCODE == i, reu_past_veg)
#     plot!(p, layer.geometry; c=ColorSchemes.turbo[i/20], opacity=0.5)
# end
# display(p)
# savefig("veg_polygon+overlap.png")


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
