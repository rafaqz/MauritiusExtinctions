using GeoJSON
using GADM
using Setfield 
using Shapefile
using RasterDataSources
using Rasters
using Rasters.LookupArrays
using Rasters: set, Between, trim, Band
using GeoInterface
using Plots
using Makie, GLMakie
using Plots: plot, plot!
using Unitful
using CategoricalArrays
includet("functions.jl")
includet("lost_land_images.jl")
years = 1638, 1773, 1835, 1872, 1935, "present"
workdir = "/home/raf/PhD/Mauritius"
datadir = joinpath(workdir, "Data")
outputdir = joinpath(datadir, "Generated")
island_keys = (; mus=:mus, reu=:reu)
borders = (
    mus=GADM.get("MUS").geom[1],
    reu=GADM.get("REU").geom[1],
    # rod=GADM.get("MUS").geom[1],
)
bbox = (
    # mus=((57.1, 57.9), (-20.6, -19.8)),
    mus=((57.1, 57.9), (-20.6, -19.949)),
    reu=((55.0, 56.0), (-22.0, -20.0)),
    # rod = (63.0, 64.0), (-20.0, -19.0),
)
tiles = getraster(SRTM; bounds=bbox.mus)
dem1 = Raster(tiles[1]; name=:DEM)
dem2 = Raster(tiles[2]; name=:DEM)
border_selectors = map(bbox) do bb
    mus=(X(Between(bb[1])), Y(Between(bb[2])), Band(1))
end
# Mauritius is right over the split in the tiles
m1 = view(dem1, border_selectors.mus...)
m2 = view(dem2, border_selectors.mus...)
mus_dem = replace_missing(trim(cat(m1, m2; dims=Y); pad=10))
# Plots.plot(mauritius_dem)
reu_tile  = getraster(SRTM; bounds=bbox.reu)[1]
reu_dem = replace_missing(read(trim(view(Raster(reu_tile), border_selectors.reu...); pad=10)))
# rod_tile  = getraster(SRTM; bounds=rod_bounds)[1]
# rod_dem = trim(view(dem3, border_selectors...); pad=10)
dems = (mus=mus_dem, reu=reu_dem)
templates = map(dems, borders) do dem, border
    mask(dem; with=boolmask(border; to=dem, shape=:polygon, boundary=:touches))
    dem
end

lc_years = 1638, 1773, 1835, 1872, 1935, "present"
lc_year_keys = map(y -> "lc_$y", lc_years)

# Flat island
# fi_selectors = X(57.64..57.69), Y(-19.9..(-19.86))
# plot(dems.mus[fi_selectors...])
# fi_pixels = length(collect(skipmissing(dems.mus[fi_selectors...])))
# step(parent(dims(dems.mus, X))) * 111
# fi_area = uconvert(u"km^2", npixels * 90.0u"m" * 90.0u"m") 

norder_dir = mkpath(joinpath(outputdir, "Norder"))
norder_rasters = replace_missing(RasterStack(norder_dir)[Band(1)])

soiltypespath = joinpath(datadir, "Norder/K factor/SoilK.shp")
soiltypenames = (Symbol.(replace.(Shapefile.Table(soiltypespath).Soil_Group, Ref(" " => "_")))...,)
soilnums = NamedTuple{soiltypenames}((1:14...,))
soilmasks = map(soilnums) do v
    norder_rasters[:soiltypes] .== v
end |> RasterStack
plot(soilmasks; c=:viridis)
plot(norder_rasters[:soiltypes])
soiltypes = norder_rasters[:soiltypes]
categorical_soil = rebuild(soiltypes; data=categorical(parent(soiltypes)), name=:categorical_soiltypes)
norder_cat_rasters = RasterStack(norder_rasters..., categorical_soil)

# lakesraster = Raster("warpedlakes.tif")[Band(1)]
# elevationraster = Raster("warpedelevation.tif")[Band(1)]
plot(dems.mus ./ maximum(skipmissing(dems.mus)))
plot!(norder_rasters[:lc_1835]; c=:viridis, opacity=0.4)
# elevationraster = replace(read(elevationraster), elevationraster[10, 10, 1] => -Inf32)

# Distance
# distance_to_coasts = map(dems, borders) do dem, border 
#     coast = boolmask(border; to=dem, shape=:line)
#     dists = nearest_distances(coast)
#     rebuild(mask(dists; with=dem); name=:to_coasts)
# end

ports = (   
     mus=(;
         major=(;
             Port_Louie=(-20.1597, 57.5012),
         ),
         minor=(;
             Port_Louie=(-20.1597, 57.5012),
             Grand_Port=(-20.3754, 57.7228),
             Grand_River_South_East=(-20.288, 57.784),
             Terre_Rouge=(-20.139, 57.499),
             Trou_dEau_Douce=(-20.244, 57.787),
             Black_River=(-20.362, 57.372),
             Poste_de_Flacq=(-20.163, 57.743),
             Grand_Baie=(-20.013, 57.584),
             Mahebourg=(-20.405, 57.709),
             Albion=(-20.218, 57.402),
             Flic_en_Flac=(-20.275, 57.371),
             Souillac=(-20.518, 57.517),
             Bel_Ombre=(-20.503, 57.399),
             Baie_du_Cap=(-20.491, 57.377),
             Poudre_de_Or=(-20.0610, 57.685),
             Grande_Gaube=(-20.009087152830663, 57.6698847150731),
             Cap=(-19.986, 57.621),
             Balaclava=(-20.0834, 57.516),
             La_Gaulette=(-20.427, 57.360),
         ),
     ),
     reu=(;
         major=(;),
         minor=(;
             Saint_Benoit=(-21.0392, 55.722),
             Le_Port=(-20.937, 55.293),
             Saint_Pierre=(-21.344, 55.481),
             Saint_Phillipe=(-21.364, 55.767),
             Saint_Paul=(-21.008, 55.270),
             Sale=(-21.269, 55.335),
             Saint_Joseph=(-21.389, 55.644),
             Saint_Denis=(-20.876, 55.446),
         ),
     ),
)

distance_to_ports = map(dems, ports) do dem, classes
    map(namedkeys(classes), classes) do key, locations
        ports = falses(dims(dem))
        for loc in locations
            selectors = Y(Contains(loc[1])), X(Contains(loc[2]))
            if DimensionalData.hasselection(ports, selectors)
                ports[selectors...] = true
            end
        end
        dists = nearest_distances(ports)
        rebuild(mask(dists; with=dem); name=Symbol("to_$(key)_ports"))
    end
end

plot(distance_to_ports.reu.minor)
plot(distance_to_ports.mus.minor)
plot(distance_to_ports.mus.major)

waterways_path = joinpath(datadir, "water.geojson")
waterways_fc = GeoJSON.read(read(waterways_path))
watermasks = map(dems) do dem 
    rebuild(boolmask(waterways_fc; to=dem); name=:waterways)
end
Plots.plot(watermasks.mus)

waterways_fc = GeoJSON.read(read(waterways_path))
watermasks = map(dems) do dem 
    rebuild(boolmask(waterways_fc; to=dem); name=:waterways)
end
Plots.plot(watermasks.mus)

distance_to_water = map(dems, watermasks) do dem, watermask
    rebuild(mask(nearest_distances(watermask); with=dem); name=:to_water)
end
plot(distance_to_water.mus)
plot(distance_to_water.reu)

# Roads
road_types = (; primary=:primary, secondary=:secondary)
highways_json = map(island_keys) do ik
    map(road_types) do rt
        json_path = joinpath(datadir, "Roads", "$(ik)_$(rt)_highways.geojson")
        GeoJSON.read(read(json_path))
    end
end
highway_masks = map(dems, highways_json) do dem, highway_type
    map(highway_type) do type
        rebuild(boolmask(type; to=dem); name=:waterways)
    end
end
distance_to_highways = map(dems, highway_masks) do dem, highways
    (; 
        primary=rebuild(mask(nearest_distances(highways.primary); with=dem); name=:to_primary_highways),
        secondary=rebuild(mask(nearest_distances(highways.primary .| highways.secondary); with=dem); name=:to_seconary_highways),
    ) 
end

Plots.plot(highway_masks.mus.secondary)
Plots.plot(distance_to_highways.mus.primary)
Plots.plot(distance_to_highways.mus.secondary)
Plots.plot(distance_to_highways.reu.primary)
Plots.plot(distance_to_highways.reu.secondary)

# Plots.plot(distance_to_water.reu; c=:seaborn_icefire_gradient, size=(1000,1000))
# Plots.plot(coasts.reu; c=:seaborn_icefire_gradient, size=(1000,1000))
# Plots.plot!(watermasks.reu)
# Plots.plot!(borders.reu; fill=nothing)
# normedelevation = 1 .- elevationraster ./ maximum(elevationraster)
# p1 = distance_to_coast .* distance_to_rivers |> plot;
# p2 = plot(elevationraster)
# p3 = plot(landuse_snapshots[5]; c=:viridis)
# plot(p1, p2, p3; layout=(1, 3))

# Slope
aspectrasters = map(d -> aspect(d, FD3Linear()), dems)
sloperasters = map(d -> slope(d, FD3Linear()), dems)
# plot(plot(sloperasters.reu), plot(aspectrasters.reu), plot(dems.reu))
# plot(plot(sloperasters.mus), plot(aspectrasters.mus), plot(dems.mus))
# plot(aspectrasters.reu)
# plot(sloperasters.mus; clims=(0, 0.5))


# Vegetation maps from "Lost Land of the Dodo"
island_rasters = map(namedkeys(island_image_classes), island_image_classes) do ik, island
    map(namedkeys(island)) do k
        read(Raster(joinpath(outputdir, "$(ik)_$(k)_cleaned.tif")))
    end
end 

# Vegetation phases as the total cleared for each time step
veg_phases = map(island_rasters, (mus=0:6, reu=0:4)) do island, phaseids
    map(phaseids) do id
        classify(island.phase, 1..id => true; others=false, missingval=missing)
    end
end

cleared_rasters = map(island_rasters) do i
    cleared = maximum(i.rem)
    (x -> x == cleared).(replace_missing(i.rem))
end

cleared_rem_plots = map(veg_phases, cleared_rasters) do vp, rem
    plot(
        plot((!).(vp[end]); c=:spring, title="Not cleared"), 
        plot((!).(rem); c=:spring, title="Remnant")
    )
end
cleared_rem_plots.reu
cleared_rem_plots.mus


# Vegetation classes
veg_classes = map(island_rasters, island_image_classes) do rasters, classes
    keys = map(x -> Symbol(replace(x, " " => "_", "-" => "_")),  classes.veg) |> values
    map((1:length(classes.veg)...,)) do id
        classify(rasters.veg, UInt8(id) => true; others=false, missingval=missing)
    end |> NamedTuple{keys} |> RasterStack
end

# Landcover
lc_dir = joinpath(datadir, "Landcover/")
# rod_tile  = RasterDataSources.zipurl(SRTM; bounds=rod_bounds)[1]
lc_shapes = (
    mus=(
        shapepath=joinpath(lc_dir, "mauritius/cla_maurice_fin.shp"),
        crspath=joinpath(lc_dir, "mauritius/cla_maurice_fin.prj"),
    ),
    reu=(
        shapepath = joinpath(lc_dir, "reunion/cla_run_2014_fin_2975.shp"),
        crspath = joinpath(lc_dir, "reunion/cla_run_2014_fin_2975.prj"),
    ),
    # rod=(
        # shapepath = joinpath(lc_dir, "rodrigues/cla_rod_fin.shp")
        # crspath = joinpath(lc_dir, "rodrigues/cla_rod_fin.prj")
    # )
)
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
lc_categories = NamedTuple{lc_names}((Int32.(0:12)...,))

# lc_rasterized_res = map(dems, lc_shapes) do dem, (; shapepath, crspath)
#     rasterize_lc(dem, shapepath, crspath; categories=lc_categories, res=90) |> read
# end
# lc_rasterized_res.mus |> plot
# lc_rasterized = map(lc_rasterized_res, dems) do A, to
#     resample(A; to) |> read
# end

# # Write to tif
# mkpath(lc_dir)
# foreach(lc_rasterized, island_keys) do lc, island
#     write(joinpath(lc_dir, "$(island)_landcover.tif"), lc)
#     rasters = Dict()
# end
# Read from tif
lc_rasterized = map(island_keys) do island
    path = joinpath(lc_dir, "$(island)_landcover.tif")
    rast =  Raster(path; name=:modern_landcover)[Band(1)]
    # Get rid of Water and NoData
    replace_missing(replace_missing(rast, lc_categories.Water), lc_categories.No_Data)
end
plot(lc_rasterized.mus)
# Masks for each land cover
lc_masks = map(lc_rasterized) do rast
    masks = map(keys(lc_categories), lc_categories) do k, v
        Raster(mask(rast .== v; with=rast, missingval=missing); name=k)
    end
    RasterStack(masks...)
end

# Count pixels of land cover categories
cover_pixels = map(lc_masks) do island
    sort([map((k, v) -> k => sum(skipmissing(v)), keys(island), values(island))...]; by=last, rev=true)
end
cover_pixels.mus
cover_pixels.reu

# Fit land cover to soil, etc
using GLM
all_layers = dems, aspectrasters, sloperasters, distance_to_water

predictor_stacks = map(soilmasks, highway_masks, lc_masks, distance_to_ports, all_layers...) do hw, lc, dp, a...
    RasterStack(lc..., dp..., a...)
end;

st = RasterStack(norder_rasters..., veg_classes.mus..., predictor_stacks.mus..., categorical_soil)
df = DataFrame(st)
plot(norder_rasters[:rainfall]/4100)
veg_classes.mus
plot!(veg_classes.mus[:semi_dry_evergreen_forest]; opacity=0.3, c=:viridis)

# model = lm(@formula(open_dry_palm_rich_woodland ~ rainfall + categorical_soiltypes + aspect + slope), df)
# r2(model)
# model = lm(@formula(semi_dry_evergreen_forest ~ DEM + categorical_soiltypes + ), df)
ag = rebuild(lc_masks.mus[:Sugarcane] .| lc_masks.mus[:Other_cropland]; name=:crops)
plot(ag; c=:viridis)
st = RasterStack(crops, soilmasks[:Latosolic_Red_Prairie_Soils], lc_masks.mus[:Barren_land], categorical_soil)
model = lm(@formula(Sugarcane ~ Latosolic_Red_Prairie_Soils), DataFrame(st))
model = lm(@formula(Barren_land ~ categorical_soiltypes), DataFrame(st))
model = lm(@formula(crops ~ categorical_soiltypes), DataFrame(st))
# model = lm(@formula(Discontinuous_urban ~ to_coasts), df)
# model = lm(@formula(Discontinuous_urban ~ to_coasts), df)
r2(model)
soilmasks

plot(soilmasks[:Latosolic_Red_Prairie_Soils]; c=:viridis)
plot(lc_masks.mus[:Discontinuous_urban]; c=:viridis)
plot(lc_masks.mus[:Discontinuous_urban] .& soilmasks[:Latosolic_Red_Prairie_Soils]; c=:viridis)
plot(crops .& soilmasks[:Latosolic_Red_Prairie_Soils]; c=:viridis)
fraction_used_urban = map(soilmasks) do A
    sum(skipmissing(lc_masks.mus[:Discontinuous_urban] .& A)) / sum(skipmissing(A))
end |> pairs |> collect
sort(fraction_used_urban; by=last, rev=true)
fraction_used_forest = map(soilmasks) do A
    sum(skipmissing((lc_masks.mus[:Forest]) .& A)) / sum(skipmissing(A))
end |> pairs |> collect
sort(fraction_used_forest; by=last, rev=true)
plot(lc_masks.mus[:Herbaceaous_vegetation] .| lc_masks.mus[:Discontinuous_urban]; c=:viridis)

fractions_by_lc = map(namedkeys(NamedTuple(lc_masks.mus))) do lckey
    s = sum(skipmissing(lc_masks.mus[lckey]))
    if s == 0
        nothing
    else
        fraction_used_ag = map(soilmasks) do A
            sum(skipmissing(lc_masks.mus[lckey] .& A)) / s
        end |> pairs |> collect |> x -> sort(x; by=x->x[2][1], rev=true)
    end
end |> pairs |> xs -> filter(x -> !isnothing(x[2]), xs)

fraction_by_soil = map(namedkeys(NamedTuple(soilmasks))) do soilkey
    map(lc_masks.mus) do A
        sum(skipmissing(soilmasks[soilkey] .& A)) / sum(skipmissing(soilmasks[soilkey]))
    end |> pairs |> collect |> x -> sort(x; by=x->x[2], rev=true)
end |> pairs

soilpixels = map(soilmasks) do A
    sum(skipmissing(A))
end |> pairs |> collect |> x -> sort(x; by=x->x[2], rev=true)

lcpixels = map(lc_masks.mus) do A
    sum(skipmissing(A))
end

sort_nt(x) = x |> pairs |> collect |> x -> sort(x; by=last, rev=true)


lc_masks.mus

lc_masks.mus[:Herbaceaous_vegetation] .| lc_masks.mus[:Discontinuous_urban]

plot(lc_masks.mus[:Herbaceaous_vegetation]; c=:viridis)
plot(lc_masks.mus[:Shrub_vegetation]; c=:viridis)
plot(lc_masks.mus[:Mangrove]; c=:viridis)

plot(lc_masks.mus[:Forest]; c=:viridis)
plot(lc_masks.mus[:Discontinuous_urban]; c=:viridis)
plot(lc_masks.mus[:Sugarcane]; c=:viridis)
plot(lc_masks.mus[:Other_cropland]; c=:viridis)
plot(lc_masks.mus[:Barren_land]; c=:viridis)

cats = lc_categories[(:Discontinuous_urban, :Forest, :Sugarcane, :Other_cropland, :Barren_land)]
plot(lc_rasterized.mus)
lc2 = clean_categories(lc_rasterized.mus;
    categories=cats, neighborhood=Moore{1}(), missingval=Int32(0),
    keep_neigborless=true, despecle=false
)
plot(lc2)
lc_masks2 = map(namedkeys(lc_categories), lc_categories) do k, v
    Raster(mask(lc2 .== v; with=lc2, missingval=missing); name=k)
end

lcpixels = map(RasterStack(lc_masks2...)) do A
    sum(skipmissing(A))
end |> sort_nt

fractions_by_lc = map(namedkeys(NamedTuple(lc_masks2))) do lckey
    s = sum(skipmissing(lc_masks2[lckey]))
    if s < 2000
        nothing
    else
        fraction_used_ag = map(soilmasks) do A
            sum(skipmissing(lc_masks2[lckey] .& A)) / s
        end |> pairs |> collect |> x -> sort(x; by=x->x[2][1], rev=true)
    end
end |> pairs |> xs -> filter(x -> !isnothing(x[2]), xs)

fraction_by_soil = map(namedkeys(NamedTuple(soilmasks))) do soilkey
    map(lc_masks2) do A
        sum(skipmissing(soilmasks[soilkey] .& A)) / sum(skipmissing(soilmasks[soilkey]))
    end |> pairs |> collect |> x -> sort(x; by=x->x[2], rev=true)
end |> pairs
