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
includet("functions.jl")
years = 1638, 1773, 1835, 1872, 1935, "present"

workdir = "/home/raf/PhD/Mauritius"
datadir = joinpath(workdir, "Data")
outputdir = joinpath(datadir, "Generated")

borders = (
    mus=GADM.get("MUS").geom[1],
    reu=GADM.get("REU").geom[1],
    # rod=GADM.get("MUS").geom[1],
)
bbox = (
    mus=((57.1, 57.9), (-20.6, -19.949)),
    reu=((55.0, 56.0), (-22.0, -20.0)),
    # rod = (63.0, 64.0), (-20.0, -19.0),
)
tiles = getraster(SRTM; bounds=bbox.mus)
dem1 = Raster(tiles[1])
dem2 = Raster(tiles[2])
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


soilraster = Raster(joinpath(outputdir, "warpedsoiltypes.tif"))[Band(1)]
Plots.plot(soilraster)
soiltypespath = joinpath(datadir, "Norder/K factor/SoilK.shp")
soiltypenames = Shapefile.Table(soiltypespath).Soil_Group
# lakesraster = Raster("warpedlakes.tif")[Band(1)]
# elevationraster = Raster("warpedelevation.tif")[Band(1)]
landuse_snapshots = map(1:6) do i
    Raster(joinpath(outputdir, "warped_landuse_snapshot_$i.tif"))[Band(1)]
end
# elevationraster = replace(read(elevationraster), elevationraster[10, 10, 1] => -Inf32)

# Distance
distance_to_coasts = map(dems, borders) do dem, border 
    coast = boolmask(border; to=dem, shape=:line)
    distance_to_coast = nearest_distances(coast)
    mask(distance_to_coast; with=dem)
end

waterways_path = joinpath(datadir, "water.geojson")
waterways_fc = GeoJSON.read(read(waterways_path))
watermasks = map(dems) do dem 
    boolmask(waterways_fc; to=dem)
end
Plots.plot(watermasks.mus)

distance_to_water = map(dems, watermasks) do dem, watermask
    mask(nearest_distances(watermask); with=dem)
end

Plots.plot(distance_to_water.mus; c=:seaborn_icefire_gradient, size=(1000,1000))
Plots.plot(coasts.reu; c=:seaborn_icefire_gradient, size=(1000,1000))
# Plots.plot!(mus_border; fill=nothing)
# normedelevation = 1 .- elevationraster ./ maximum(elevationraster)
# p1 = distance_to_coast .* distance_to_rivers |> plot;
# p2 = plot(elevationraster)
# p3 = plot(landuse_snapshots[5]; c=:viridis)
# plot(p1, p2, p3; layout=(1, 3))

# Slope
aspectrasters = map(d -> aspect(d, FD3Linear()), dems)
sloperasters = map(d -> slope(d, FD3Linear()), dems)
Plots.plot(Plots.plot(sloperasters.reu), Plots.plot(aspectrasters.reu), Plots.plot(dems.reu))
Plots.plot(Plots.plot(sloperasters.mus), Plots.plot(aspectrasters.mus), Plots.plot(dems.mus))
Plots.plot(sloperasters.reu; clims=(0,0.5))
Plots.plot(aspectrasters.reu)
sloperasters.reu
aspectrasters.reu

ps = map(r -> Plots.plot(r; c=:terrain, size=(1000, 1000), clims=(0, 1.0)), sloperasters)
Plots.plot(ps...)

lc_lookup = Dict(lc_categories .=> 0:12)

lc_path = joinpath(datadir, "Landcover/")
# rod_tile  = RasterDataSources.zipurl(SRTM; bounds=rod_bounds)[1]
lc_shapes = (
    mus=(
        shapepath=joinpath(lc_path, "mauritius/cla_maurice_fin.shp"),
        crspath=joinpath(lc_path, "mauritius/cla_maurice_fin.prj"),
    ),
    reu=(
        shapepath = joinpath(lc_path, "reunion/cla_run_2014_fin_2975.shp"),
        crspath = joinpath(lc_path, "reunion/cla_run_2014_fin_2975.prj"),
    ),
    # rod=(
        # shapepath = joinpath(lc_path, "rodrigues/cla_rod_fin.shp")
        # crspath = joinpath(lc_path, "rodrigues/cla_rod_fin.prj")
    # )
)

lc_rasterized = map(dems, lc_shapes) do dem, (; shapepath, crspath)
    rasterize_lc(dem, shapepath, crspath)
end

lc_dir = joinpath(outputdir, "Landcover")
mkpath(lc_dir)
island_keys = NamedTuple{keys(lc_shapes)}(keys(lc_shapes))
# Masks for each land cover
foreach(lc_rasterized, island_keys) do lc, island
    write(joinpath(lc_dir, "$(island)_landcover.tif"), lc)
    rasters = Dict()
    for (k, v) in lc_lookup
        mask = Raster(UInt8.(lc .== v); missingval=nothing)
        layername = replace(k, " " => "_")
        write(joinpath(lc_dir, "$(island)_$(layername).tif"), mask)
    end
end

lc_rasterized = map(island_keys) do island
    return Raster(joinpath(lc_dir, "$(island)_landcover.tif"))
end

lc_masks = map(island_keys) do island
    masks=Dict()
    for (k, v) in lc_lookup
        layername = replace(k, " " => "_")
        masks[k] = Bool.(Raster(joinpath(lc_dir, "$(island)_$(layername).tif")))
    end
    return masks
end

plot_lc_makie(lc_rasters.reu.landcover)
Makie.heatmap(parent(read(lc_rasterized.reu
                          [Band(1)])))
Plots.heatmap(read(lc_masks.reu[Band(1)]))
plot_lc(lc_rasters.reu.landcover)
Plots.plot(lc_rasters.reu.masks["Sugarcane"]; c=:spring)
Plots.plot(lc_rasters.reu.masks["Forest"]; c=:spring)
Plots.plot(lc_rasters.reu.masks["Shrub vegetation"]; c=:spring)
Plots.plot(lc_rasters.reu.masks["Discontinuous urban"]; c=:spring)
Plots.plot(lc_rasters.reu.masks["Other cropland"]; c=:spring)
keys(lc_rasters.mus.masks)

cover_pixels = map(lc_rasters) do island
    Dict(map((k, v) -> k => sum(v), keys(island.masks), values(island.masks)))
end
