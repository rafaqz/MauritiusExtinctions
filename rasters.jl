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
soiltypes = norder_stack[:soiltypes]
categorical_soil = rebuild(soiltypes; data=categorical(parent(soiltypes)), name=:categorical_soiltypes)

# lakesraster = Raster("warpedlakes.tif")[Band(1)]
# elevationraster = Raster("warpedelevation.tif")[Band(1)]
plot(dems.mus ./ maximum(skipmissing(dems.mus)))
plot!(norder_stack[:lc_1835]; c=:viridis, opacity=0.4)
# elevationraster = replace(read(elevationraster), elevationraster[10, 10, 1] => -Inf32)

distance_stack = map(island_keys) do i
    read(RasterStack(joinpath(distancedir, string(i)))[Band(1)])
end
distance_stack.mus

# Slope
slope_stack = map(dems) do dem
    slopeaspect(dem, FD3Linear()) 
end
plot(slope_stack.reu)

# Vegetation maps from "Lost Land of the Dodo"
lostland_stack = map(namedkeys(lostland_image_classes), lostland_image_classes) do i, rasters
    read(RasterStack(joinpath(outputdir, "LostLand", string(i))))
end 
plot(lostland_stack.reu)

# Vegetation classes
veg_classes = map(lostland_stack, lostland_image_classes) do stack, classes
    map(stack, classes) do r, c
        keys = map(x -> Symbol(replace(x, " " => "_", "-" => "_")), c) |> values
        map((1:length(c)...,)) do id
            classify(r, UInt8(id) => true; others=false, missingval=missing)
        end |> NamedTuple{keys} |> RasterStack
    end
end

plot(veg_classes.mus.veg; c=:viridis)

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
