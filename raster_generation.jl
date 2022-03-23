using GeoJSON
includet("raster_common.jl")

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

# To ports
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

# To water
waterways_path = joinpath(datadir, "water.geojson")
waterways_fc = GeoJSON.read(read(waterways_path))
watermasks = map(dems) do dem 
    rebuild(boolmask(waterways_fc; to=dem); name=:waterways)
end

distance_to_water = map(island_keys, dems, watermasks) do i, dem, watermask
    rast = mask(nearest_distances(watermask); with=dem)
    write(joinpath(distancedir, string(i), "to_water.tif"), rast)
    rast
end

# To roads
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
distance_to_roads = map(island_keys, dems, highway_masks) do i, dem, highways
    rast = mask(nearest_distances(highways.primary); with=dem)
    write(joinpath(distancedir, string(i), "to_primary_roads.tif"), rast)
    rast = mask(nearest_distances(highways.primary .| highways.secondary); with=dem)
    write(joinpath(distancedir, string(i), "to_secondary_roads.tif"), rast)
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
