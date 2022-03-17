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
    mus=((57.1, 57.9), (-20.6, -19.8)),
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
plot(dems.reu)

# Flat island
fi_selectors = X(57.64..57.69), Y(-19.9..(-19.86))
plot(dems.mus[fi_selectors...])
fi_pixels = length(collect(skipmissing(dems.mus[fi_selectors...])))
step(parent(dims(dems.mus, X))) * 111
fi_area = uconvert(u"km^2", npixels * 90.0u"m" * 90.0u"m") 


ports = (
     mus=(
         Port_Louie=((-20.1597, ations57.5012), 100),
         Grand_Port=((-20.3754, 57.7228), 30),
         Grand_River_South_East=((-20.288, 57.784), 10),
         # Black_River=((-20.362, 57.372), 10),
         # Flacq=((-20.243, 57.7850), 5),
     ),
     reu=(
     ),
)


soilraster = Raster(joinpath(outputdir, "warpedsoiltypes.tif"))[Band(1)]
Plots.plot(soilraster)
soiltypespath = joinpath(datadir, "Norder/K factor/SoilK.shp")
soiltypenames = Shapefile.Table(soiltypespath).Soil_Group
# lakesraster = Raster("warpedlakes.tif")[Band(1)]
# elevationraster = Raster("warpedelevation.tif")[Band(1)]
landuse_snapshots = map(1:6) do i
    Raster(joinpath(outputdir, "warped_landuse_snapshot_$i.tif"))[Band(1)]
end
plot(dems.mus ./ maximum(skipmissing(dems.mus)))
plot!(landuse_snapshots[2]; c=:viridis, opacity=0.4)
# elevationraster = replace(read(elevationraster), elevationraster[10, 10, 1] => -Inf32)

# Distance
distance_to_coasts = map(dems, borders) do dem, border 
    coast = boolmask(border; to=dem, shape=:line)
    dists = nearest_distances(coast)
    mask(dists; with=dem)
end

distance_to_ports = map(dems, ports) do dem, locations
    ports = falses(dims(dem))
    for loc in locations
        selectors = Y(Contains(loc[1][1])), X(Contains(loc[1][2]))
        if DimensionalData.hasselection(ports, selectors)
            ports[selectors...] = true
        end
    end
    dists = nearest_distances(ports)
    mask(dists; with=dem)
end
plot(distance_to_ports.mus)
plot(distance_to_ports.reu)

waterways_path = joinpath(datadir, "water.geojson")
waterways_fc = GeoJSON.read(read(waterways_path))
watermasks = map(dems) do dem 
    boolmask(waterways_fc; to=dem)
end
Plots.plot(watermasks.mus)

distance_to_water = map(dems, watermasks) do dem, watermask
    mask(nearest_distances(watermask); with=dem)
end

Plots.plot(distance_to_water.reu; c=:seaborn_icefire_gradient, size=(1000,1000))
Plots.plot(coasts.reu; c=:seaborn_icefire_gradient, size=(1000,1000))
Plots.plot!(watermasks.reu)
Plots.plot!(borders.reu; fill=nothing)
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
Plots.plot(aspectrasters.reu)
Plots.plot(sloperasters.reu; clims=(0, 0.2))

# Landcover
lc_lookup = Dict(lc_categories .=> 0:12)

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

# lc_rasterized = map(dems, lc_shapes) do dem, (; shapepath, crspath)
    # rasterize_lc(dem, shapepath, crspath)
# end

# mkpath(lc_dir)
island_keys = NamedTuple{keys(lc_shapes)}(keys(lc_shapes))
# # Masks for each land cover
# foreach(lc_rasterized, island_keys) do lc, island
#     write(joinpath(lc_dir, "$(island)_landcover.tif"), lc)
#     rasters = Dict()
#     for (k, v) in lc_lookup
#         mask = Raster(UInt8.(lc .== v); missingval=nothing)
#         layername = replace(k, " " => "_")
#         write(joinpath(lc_dir, "$(island)_$(layername).tif"), mask)
#     end
# end

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

Plots.plot(lc_masks.reu["Sugarcane"]; c=:spring)
Plots.plot(lc_masks.reu["Forest"]; c=:spring)
Plots.plot(lc_masks.reu["Shrub vegetation"]; c=:spring)
Plots.plot(lc_masks.reu["Discontinuous urban"]; c=:spring)
Plots.plot(lc_masks.reu["Other cropland"]; c=:spring)
Plots.plot(lc_masks.mus["Sugarcane"]; c=:spring)
Plots.plot(lc_masks.mus["Forest"]; c=:spring)
Plots.plot(lc_masks.mus["Shrub vegetation"]; c=:spring)
Plots.plot(lc_masks.mus["Discontinuous urban"]; c=:spring)
Plots.plot(lc_masks.mus["Other cropland"]; c=:spring)
keys(lc_rasters.mus.masks)

# Count pixels of land cover categories
cover_pixels = map(lc_masks) do island
    sort(map((k, v) -> k => sum(v), keys(island), values(island)); by=last, rev=true)
end
cover_pixels.mus
cover_pixels.reu
