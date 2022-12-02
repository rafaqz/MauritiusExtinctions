includet("raster_common.jl")
includet("ports.jl")

# Distances

# To coasts
distance_to_coasts = map(dems, borders) do dem, border 
    coast = boolmask(border; to=dem, shape=:line)
    dists = nearest_distances(coast)
    rebuild(mask(dists; with=dem); name=:to_coasts)
end
foreach(island_keys, distance_to_coasts) do i, raster
    mkpath(joinpath(distancedir, string(i)))
    write(joinpath(distancedir, string(i), "to_coast.tif"), raster)
end

foreach(island_keys, dems, ports) do i, dem, classes
    foreach(namedkeys(classes), classes) do k, locations
        ports = falses(dims(dem))
        for loc in locations
            selectors = Y(Contains(loc[1])), X(Contains(loc[2]))
            if DimensionalData.hasselection(ports, selectors)
                ports[selectors...] = true
            end
        end
        dists = nearest_distances(ports)
        rast = rebuild(mask(dists; with=dem); name=Symbol("to_$(k)_ports"))
        dir = joinpath(distancedir, string(i))
        write(joinpath(dir, "to_$(k)_ports.tif"), rast)
    end
end

distance_to_water = map(island_keys, dems, watermasks) do i, dem, watermask
    rast = mask(nearest_distances(watermask); with=dem)
    write(joinpath(distancedir, string(i), "to_water.tif"), rast)
    rast
end

# Landcover
lc_dir = joinpath(datadir, "Landcover/")
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
lc_categories = NamedTuple{lc_names}((Int32.(0:12)...,))

lc_rasterized_res = map(dems, lc_shapes) do dem, (; shapepath, crspath)
    rasterize_lc(dem, shapepath, crspath; categories=lc_categories, res=90) |> read
end
lc_rasterized_res.mus |> plot
lc_rasterized = map(lc_rasterized_res, dems) do A, to
    resample(A; to) |> read
end

# Write to tif
mkpath(lc_dir)
foreach(lc_rasterized, island_keys) do lc, island
    write(joinpath(lc_dir, "$(island)_landcover.tif"), lc)
    rasters = Dict()
end


ag = Raster("/home/raf/PhD/Mauritius/Data/Deforestation-Mauritius/rasters/agr_suitab/agr_suitab_h/w001001.adf"; crs=EPSG(3337))
plot(ag)
resample(ag; to=dems.mus) |> plot
