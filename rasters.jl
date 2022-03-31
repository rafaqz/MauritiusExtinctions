includet("raster_common.jl")

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

distance_stacks = map(island_keys) do i
    read(RasterStack(joinpath(distancedir, string(i)))[Band(1)])
end
distance_stacks.mus

# Slope
slope_stacks = map(dems) do dem
    slopeaspect(dem, FD3Linear()) 
end
plot(slope_stacks.reu)

# Vegetation maps from "Lost Land of the Dodo"
lostland_stacks = map(namedkeys(lostland_image_classes), lostland_image_classes) do i, rasters
    read(RasterStack(joinpath(outputdir, "LostLand", string(i))))
end 
plot(lostland_stacks.reu)

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
