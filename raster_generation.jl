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
# In both cases we have a shift in the
# location of the major port.
ports = (   
    mus=(;
        major1=(;
            Grand_Port=(-20.3754, 57.7228),
        ),
        major2=(;
            Port_Louie=(-20.1597, 57.5012),
        ),
        minor=(;
            Albion=(-20.218, 57.402),
            Baie_du_Cap=(-20.491, 57.377),
            Bel_Ombre=(-20.503, 57.399),
            Black_River=(-20.362, 57.372),
            Cap=(-19.986, 57.621),
            Flic_en_Flac=(-20.275, 57.371),
            Grand_Baie=(-20.013, 57.584),
            Grand_Port=(-20.3754, 57.7228),
            Grand_River_South_East=(-20.288, 57.784),
            Grande_Gaube=(-20.009, 57.670),
            La_Gaulette=(-20.427, 57.360),
            Mahebourg=(-20.405, 57.709),
            Port_Louie=(-20.1597, 57.5012),
            Poste_de_Flacq=(-20.163, 57.743),
            Poudre_de_Or=(-20.0610, 57.685),
            Souillac=(-20.518, 57.517),
            Terre_Rouge=(-20.139, 57.499),
            Trou_dEau_Douce=(-20.244, 57.787),
            Turtle_Bay=(-20.0834, 57.516),
        ),
    ),
    reu=(;
        major1=(;
            Saint_Denis=(-20.876, 55.446),
        ),
        major2=(;
            Le_Port=(-20.937, 55.293),
        ),
        minor=(;
            Le_Port=(-20.937, 55.293),
            Saint_Paul=(-21.008, 55.270),
            Saint_Benoit=(-21.0392, 55.722),
            Saint_Denis=(-20.876, 55.446),
            Saint_Joseph=(-21.389, 55.644),
            Saint_Pierre=(-21.344, 55.481),
            Saint_Phillipe=(-21.364, 55.767),
            Sale=(-21.269, 55.335),
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
    json_path = joinpath(datadir, "Roads", "$(ik)_highways.geojson")
    GeoJSON.read(read(json_path))
end

function namelike(names, feature)
    haskey(feature.properties, "name") || return false
    return any(occursin(feature.properties["name"]), names)
end

# 1764 - Belin 
# "Port Louis – The South Motorway"
# "Brabant Street"
# "Royal Road Grand River"
# "Royal Road Beau Bassin"
# "Royal Road Rose Hill"
# "Royal Road Belle Rose"
# "Royal Road Phoenix"
# "Royal Road Curepipe"
# "Royal Road Eau Coulee"
# "A10"
# "Royal Road Castel"
# "Royal Road Saint Paul"
# "Route Royale Phoenix – Plaisance"
# "Brisee Verdiere - Saint Julien - Constance Road"
# "Brisee Verdiere-Saint Julien-Constance Road"
# "Old Flacq Road"
# "Plaissance Ferney Road"
# "Route Royale Flacq – Mahebourg"
# "Route Royale Plaine Magnien – Mahebourg"
# "Royal Road Forest Side"
# "Grand River South East Road"
# "Port Louis - Central Flacq Road"
# "Pamplemousses Road"
# "Royal Street"
# "Lord Kitchener Street"
# "John Kennedy Street"
# "Old Moka Road"
# "Pont Fer Flyover"
# # "Moka Camp de Masque Flacq Road"
# # 1812 Melbert
# "A9" 
# "B9" 
# "Royal Road La Flora"
# "Savanne Road"
# "Royal Road Britannia"
# "Royal Road Tyack"
# "Royal Road Rivière des Anguilles"
# "Royal Road Saint Aubin"
# "Royal Road Souillac"
# "Coastal Road Souillac"
# "Coastal Road Surinam"
# "Coastal Road Black River – Savanne"
# "Coastal Road Rivière des Galets"
# "Coastal Road Saint Martin"
# "Coastal Road Baie du Cap"
# "Coastal Road Le Morne"
# "Coastal Road La Gaulette"
# "Coastal Road Case Noyale"
# "Coastal Road Small Black River"
# "Coastal Road Black River"
# "Royal Road Black River"
# "Royal Road Tamarin"
# "Royal Road Bambous"
# "Royal Road Canot"
# "Royal Road Gros Cailloux"
# "Royal Road Richelieu"
# "Royal Road Petite Rivière"
# "La Barraque Road"
# "La Baraque Road"
# "Riche Bois Road"
old_road_names = [
 "A1"
 "A10"
 "A12"
 "A13"
 "A15"
 "A16"
 "A2"
 "A3"
 "A4"
 "A5"
 "A6"
 "A7"
 "A8"
 "A9"
 "B10"
 "B127" # Near M1
 "B128" # Near M1
 "B15"
 "B2"
 "B23"
 "B24"
 "B26"
 "B27"
 "B28"
 "B3"
 "B48"
 "B5"
 "B52"
 "B55"
 "B61"
 "B7"
 "B70"
 "B8"
 "B83"
 "B84"
 "B86"
 "B89"
 "B9"
 "B90"
 "B94"
 "M2"
 "Vacoas la Marie Road"
]

function roadin(feature, roadlist)
    props = feature.properties
    haskey(props, "highway") || return false
    props["highway"] in ("primary", "secondary", "tertiary", "unclassified", "residential") || return false
    return get(props, "ref", "") in roadlist ||
        get(props, "name", "") in roadlist ||
        get(props, "@id", "") in roadlist
end

road_names_1725 = [
   "A1"
   # Petite Riviere/Albion
   "Royal Road Petite Rivière"
   "Albion Road"
   # "A3"
   "Royal Road Richelieu"
   # Pamplemousses
   "A2"
   # Flacq
   "way/22956430" # "B23"
   # Grand River South East
   "B28"
   "Grand River South East Road"
   # Savanne
   "Coastal Road Souillac"
   "Coastal Road Surinam"
   "Coastal Road Rivière des Galets"
   "Coastal Road Black River – Savanne" # Maybe not the Western end?
]

road_names_1764 = [road_names_1725...
   "A10"
   "Pont Fer Flyover"
   "A6"
   # Port Louie
   "B143"
   # Pamplemousses
   "M2"
   "B39"
   "B20"
   "B119"
   "Trou aux Biches Road"
   "Terre Rouge Triolet vGrand Baie Road"
   "Terre Rouge Triolet Grand Baie Road"
   "B97"
   # Grand Port
   "B28"
   # "A12"
   "A15"
   # Moka
    "A7"
   "B54"
   "A17"
   "B27"
   # Flaq
   "B56"
   "B60"
   "B23"
   "B17"
   # 
]

road_names_1813 = [road_names_1764...
    # North
    "B13"
    "B15"
    # Flacq
    "B27" # Three islots
    # Moka
    "B46"
    "B49"
    "B93"
    # Grand Port
    "B8"
    "B79"
    "B79"
    "B84"
    # Savanna
    "B102"
    "B89"
    "B9"
    "A3"
]

# old_road_names = road_names_1813
old_road_names = road_names_1764
highways = highways_json.mus.features
selected_roads = [f for f in highways if roadin(f, old_road_names)]
selected_roads
# plot(dems.mus)
plot(deforestation.mus[:by_1807]; c=:viridis)
plot!(selected_roads)

plot(distance_stacks.mus[:to_minor_ports])

for r in old_road_names
    r in all_road_names || println(r)
end
all_road_names = sort([f.properties["name"] for f in highways if haskey(f.properties, "name")])
filter(n -> occursin("Triolet", n), all_road_names)

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
