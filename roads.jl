includet("raster_common.jl")
using Dates
using DimensionalData.LookupArrays

# To roads
function namelike(names, feature)
    haskey(feature.properties, "name") || return false
    return any(occursin(feature.properties["name"]), names)
end
function roadin(feature, roadlist)
    props = feature.properties
    haskey(props, "highway") || return false
    props["highway"] in ("trunk", "primary", "secondary", "tertiary", "unclassified", "track", "residential") || return false
    return get(props, "ref", "") in roadlist ||
        get(props, "name", "") in roadlist ||
        get(props, "@id", "") in roadlist
end
select_ways(ways, way_names) = [f for f in ways if roadin(f, way_names)]
function plot_ways(ways, way_names, dem, title="")
    selected_ways = select_ways(ways, way_names)
    plot(dem ./ maximum(skipmissing(dem)), legend=false)
    # plot((dems.mus ./ maximum(skipmissing(dems.mus)))[X=57.4..58.6, Y=(-20.35)..(-20.15)])
    # plot!(watermasks.mus; color=:blue, legend=false) 
    # plot!(deforestation.mus[:by_1835]; c=:viridis, opacity=0.4) 
    # plot!(selected_ways; color=:red, title, legend=false)
    plot!(selected_ways; title, legend=false)
end

ways_json = map(island_keys) do ik
    json_path = joinpath(datadir, "Roads", "$(ik)_ways.geojson")
    GeoJSON.read(read(json_path))
end

waterways_path = joinpath(datadir, "water.geojson")
waterways_fc = GeoJSON.read(read(waterways_path))
watermasks = map(dems) do dem 
    rebuild(boolmask(waterways_fc; to=dem); name=:waterways)
end

mus_way_names_1725 = [
    # Port Louis
    "Royal Road Grand River"
    "Brabant Street"
    "Lord Kitchener Street"
    # Petite Riviere/Albion
    "Royal Road Petite Rivière"
    # "Albion Road"
    # "A3"
    "Royal Road Richelieu"
    # Pamplemousses
    "A2" # The middle hump has changed to go around rather than through the Flaq rivers.
    # Flacq
    "B117"
    "way/22956430" # "B23"
    # Grand River South East
    "Grand River South East Road"
    # "Route Royale Flacq – Mahebourg"
    "way/22956445"
    "way/255665959"
    "way/203228490"
    # Savanne
    "Coastal Road Souillac"
    "Coastal Road Surinam"
    "Coastal Road Rivière des Galets"
]
mus_way_names_1764 = [mus_way_names_1725...
    # Port Louie
    "Military Road"
    # "Queen Street" # Has dups
    "David Street"
    "B32"
    "B77"
    "B135"
    "B136"
    "B137"
    "B139"
    "B143"
    "Seychelles Street"
    # "Mahatma Gandhi Street" # Has dups
    "Pouce Street"
    # "Sir Edgar Laurent Street" # Has dups
    "way/22821395"
    "way/372815956"
    "way/499140220"
    "way/509015472"
    "way/172973974"
    "way/172973971"
    "way/178752989"
    # Plains de Willhelm
    "A1"
    # Pamplemousses
    "B18"
    "B20"
    "B39"
    # "A4" # But only the lower part
    "Terre Rouge Triolet vGrand Baie Road"
    "Trou aux Biches Road"
    "way/97536024"
    # Rivier du Rempart
    "B21"; "way/385138585"; "way/464693923"; "way/385138589";
    "B160"
    "B161"
    # A6 
    # Grand Port
    # Nearest part of B28
    "way/111531369"
    "way/632723514"
    # "way/607624541" # Has too much going east 
    "way/172366496"
    # Old North-South road over the mountain
    "way/444653257"
    "way/500756668"
    "way/761226795"
    # "way/444653259" # Newer?
    # Moka
    # Old Moka Grand Port Road
    "way/515144171" # A7 to Militair
    "way/422974133"
    "way/934085659"
    "way/422974132"
    # "way/269927546"
    "Dubreuil to Mt Blanche track"
    # "way/422974130" # Has extra road
    "Vullemin Link Road"
    # "way/408737644" # To far up the valley
    "way/408737638"
    "way/408737648"
    "way/422974124"
    "B47"
    # "B54" dups
    "A17"
    # Flaq
    "way/22956438" # Part of A7
    # "B17"
    "B23"
    "B26"
    "B60"
]
mus_way_names_1813 = [mus_way_names_1764...
    # North
    "B13"
    # Pamplemousses
    "M2"
    "B119"
    "Trou aux Biches Road"
    "Terre Rouge Triolet vGrand Baie Road"
    "Terre Rouge Triolet Grand Baie Road"
    "B97"
    # Flacq
    "B27" # Three islots
    "B56"
    "way/982769094"
    "way/750288522"
    "Montaigu Road"
    "Petit Bois Road"
    "Bon Accueil Branch Road"
    "way/173830685" # Camp Thorel
    "way/297545245"
    "Nouvelle Découverte Road"
    "Unité Road"
    # Sebastopol has unnamed roads
    "Petit Paquet Road"
    "way/256693691"; "Pellegrin Road"; "Sebastopol Road"; "way/192695526"; "way/23015209";
    "Sans Souci Road"
    # Moka
    # M1
    "way/158659416"; "way/688991901"; "way/1000761137"; "way/688991907"; "way/22955614";
    "A7"
    "B24"
    "B27"
    "B46"
    "B49" # But not the end near the dam
    "B52"
    "B96"
    "B93"
    "B169"
    "Dubreuil to Mt Blanche track"
    "Quartier Militaire Road"
    "Cote d'Or Road"
    # Plains de Willhelm
    "A10"
    "Mare Longue track"
    "way/510324249"
    "way/511272692"
    "way/510324248"
    "way/71499042"
    "way/71499043"
    "way/1005885837"
    "way/1018759143"
    "way/260320939"
    "way/260320939"
    "way/71499042"
    "way/186185676"
    "way/85416827"
    "Pont Fer Flyover"
    "way/510027720"
    "way/997358228"
    "way/510027720"
    "way/22955617"
    "way/22955618"
    "A6"
    # Grand Port
    "A9"
    "B8"
    "B28"
    "B82"
    "B84"
    "A12"
    "A15"
    "way/234589175"
    "way/22980286"
    "way/22980286"
    "Shivala Road"
    # East-West mid road
    "way/525674976"
    "way/587611655"
    "way/1007800223"
    "way/473933406"
    "way/587611654"
    "way/473932236"
    # Above Plain Margien
    "way/111531365"
    "B7"
    "B83"
    # Joins south St Hubert, but is actually 2 blocks above the original
    "way/954499321"
    "way/348699426" 
    # Connection north
    "way/464277319" 
    "way/731189680" # Not the actual road - east of it. The original is gone
    "way/731189678"
    "way/731189677"
    "way/731183127"
    "way/731183127"
    "way/731183128"
    "way/683712368"
    "way/683712365"
    "way/683712366"
    # Savanne
    "B89"
    "B9"
    "A3"
    # Black River
    "B2"
    "B31" # Except it curves to far north now, old southern end is abandonned
    "B34"
    "B91"
    "Yemen Road"
    "Black River Yemen Road"
    "B114"
    "Les Salines Pilot Road"
    "Les Salines Pilote Road"
    "way/116709096"
    "way/729188446"
    "way/402006019"
    "way/729191999"
    "way/270046786"
]
mus_way_names_1880 = [mus_way_names_1813...
    # Port Louie
    "B19"
    "B34"
    "way/178844181"
    "way/858702316"
    "way/499957228"
    "way/388705835"
    "way/178844182"
    "way/499957228"
    "way/87081541"
    # Rivier du Rempart
    "A5"
    "B11"
    "B14"
    "B15" # But not the hump
    "B16"
    "B17"
    "B44"
    "B45"
    "B171"
    "Union Road"
    "B100"
    "Forbach Road"
    "Fond du Sach Forbach Road"
    "Chemin Hollandais Coast Road"
    # Pamplemousses
    "Arsenal Balaclava Road"
    "B29"
    "B35"
    "B38"
    # Flacq
    "B22"
    "B59"
    "B55"
    "B61"
    "B62"
    "B99"
    # Moka
    # Grand Port
    "A12"
    "Deux Bras – Cent Gaulettes Road"
    "B81"
    "Mon Desert Road"
    "way/379773482" # North end
    "B88" # But just the start
    "B90"
    "Saint Avold Street"; "way/22980335"; "way/22980328"; "way/22953306";
    "way/22953254"
    # Savanne
    "B102"
    "way/260900945"
    "way/260900951"
    "way/260900950"
    "way/260904675"
    "Bon Courage"
    "way/1022927355"
    "B104"
    # Plains de Willhelm
    "A8"
    "B3"
    "Vacoas la Marie Road"
    "Beau Plateau Branch Road"
    "way/211462413"
    "way/859758776"
    "way/252646178"
    "B43"
    "B4"
    "B11"
    "B17"
    "B37"
    "B68"
    "B70"
    "B130"
    "B167"
    "way/71499043"
    "Valentina – Bagatelle Link Road (Northbound)"; "Valentina – Bagatelle Link Road"; "way/407080838"; "way/407080851";
    "Pandit  J.Nehru Road"
    # Black River
    "B103"
    "Anna Branch Road"
    "way/499161502"
    "way/319123850"
    "Magenta Road"
    "way/407269865"
    "way/407506295"
    "way/407506298"
    "way/407506296"
]


reu_way_names_1752 = [
    # Center
    "N 3"
]
reu_way_names_1752 = [
    # Center
    "N 3"
]
reu_way_names_1804 = [
    # "Chemin Agenor"
    # South
    "N 2"
    # Center
    "N 3"
    # Cilaos
    "N 5"
    # Salazie
    "D 48"
    # North
    "N 1"
    "N 1A"
    "N 1a"
    "N 1C"
    "N 1E"
    # "Route du Littoral" # Broken??
    # "N 1"
    # "N 4"
    # "N 6"
    # "N 7"
    # "N 8"
]

mus_way_names = mus_way_names_1725
# mus_way_names = mus_way_names_1764
# mus_way_names = mus_way_names_1813
# mus_way_names = mus_way_names_1880
reu_way_names = reu_way_names_1804
ways = map(x -> x.features, ways_json)
way_names = (mus=mus_way_names, reu=reu_way_names)
plot(dems.reu)
plot!(select_ways(ways.mus, way_names.mus))
# plot_ways(ways.reu, way_names.reu, dems.reu)

ag = Raster("/home/raf/PhD/Mauritius/Data/Deforestation-Mauritius/rasters/agr_suitab/agr_suitab_h/w001001.adf"; crs=EPSG(3337))
plot(ag)
resample(ag; to=dems.mus) |> plot

mus_way_names_sequence = way_names_1725, way_names_1764, way_names_1813, way_names_1880
mus_way_years = (1725, 1764, 1813, 1880)
# plot(plot_ways.(Ref(ways), way_names_sequence, string.(way_years))...)
# savefig("roads_by_year.png")


all_way_names = map(ways) do w
    sort([f.properties["name"] for f in w if haskey(f.properties, "name")])
    # sort([f.properties["ref"] for f in w if haskey(f.properties, "ref")])
end
# names = filter(n -> occursin("Littoral", n), all_way_names.reu)

function roadin(feature, roadlist)
    props = feature.properties
    haskey(props, "highway") || return false
    props["highway"] in ("trunk", "primary", "secondary", "tertiary", "unclassified", "track", "residential") || return false
    return get(props, "ref", "") in roadlist ||
        get(props, "name", "") in roadlist# ||
        get(props, "@id", "") in roadlist
end
select_ways(ways.reu, names)



function rasterize_ways(ways, way_names; to)
    selected_ways = select_ways(ways, way_names)
    rasterize(selected_ways; to, fill=true) 
end

rasterize_ways(ways, way_names_1725; to=dems.mus)

selected_ways = select_ways(ways, way_names_1880)
boolmask(selected_ways; to=dems.mus, shape=:line)


way_masks = map(way_names_sequence, road_years) do way_names, year
    selected_ways = select_ways(ways, way_names)
    Raster(boolmask(selected_ways; to=dems.mus, shape=:line); name = Symbol("ways_$year"))
end |> RasterStack
way_distances = map(nearest_distances, way_masks)
plot(way_masks)
plot(way_distances; clims=(0, 320))
way_distances

to_ways = Raster(way_distances; dims=Ti(Date.([road_years...]); sampling=Intervals(Center())))
plot(to_ways)

# way_masks = map(dems, ways_json) do dem, way_type
#     map(way_type) do type
#         rebuild(boolmask(type; to=dem); name=:waterways)
#     end
# end



plot(distance_stacks.mus[:to_minor_ports])

for r in old_road_names
    r in all_road_names || println(r)
end
all_road_names = sort([f.properties["name"] for f in ways if haskey(f.properties, "name")])
filter(n -> occursin("Valentina", n), all_road_names)
filter(f -> occursin("way/510027720", f.properties["@id"]), ways)


distance_to_roads = map(island_keys, dems, way_masks) do i, dem, ways
    rast = mask(nearest_distances(ways.primary); with=dem)
    write(joinpath(distancedir, string(i), "to_primary_roads.tif"), rast)
    rast = mask(nearest_distances(ways.primary .| ways.secondary); with=dem)
    write(joinpath(distancedir, string(i), "to_secondary_roads.tif"), rast)
end
