using JSON3, MapRasterization, GeoInterface

function get_map_files()
    selected_dir = "/home/raf/PhD/Mascarenes/Data/Selected"
    file_details = (
        mus=(;
        atlas_1992_agriculture = (filename="Mauritius/Undigitised/atlas_1992_agriculture.jpg", poly=1, layers=(;
            urban="urban",
            forest="forest",
            forestry="forestry",
            cleared=(|, "cane", "forage", "tea", "market_gardens", "forestry", "urban", "cane_or_market_gardens", "cane_or_tea", "tea_or_market_gardens", "cane_or_fruit", "pasture", "lakes"),
        )),
        atlas_1992_land_use = (filename="Mauritius/Undigitised/atlas_1992_land_use.jpg", poly=1, layers=(;
            urban=(|, "urban", "other_state_land_urban"),
            forest=(|, "forest", "mountain_reserves", "tea_development_forest", "private_forest_or_wasteland"),
            agriculture=(|, "small_properties", "medium_properties", "large_properties", "rose_bell"),
            cleared=(|, "urban", "other_state_land_urban", "other_state_land", "small_properties", "medium_properties", "large_properties", "rose_bell"),
        )),
        vegetation = (filename="Mauritius/Undigitised/page33_mauritius_vegetation_colored.png", poly=1, layers=(;)),
        atlas_dutch_period = (filename="Mauritius/Undigitised/atlas_dutch_period.jpg", poly=1, layers=(;
            uncleared="undisturbed",
            ebony_harvest="ebony_harvest",
            cleared="cleared",
        )),
        atlas_18C_land_use = (filename="Mauritius/Undigitised/atlas_18C_land_use.jpg", poly=1, layers=(;
            cleared=(
                cleared_1772=(|, "urban_1763", "cleared_1772", "abandonned_1810"),
                cleared_1810=(|, "urban_1763", "cleared_1772", "urban_1810", "cleared_1810"),
            ),
            urban = (urban_1763="urban_1763", urban_1810=(|, "urban_1763", "urban_1810")),
            abandonned=(abandonned_1810="abandonned_1810",),
            # uncleared=(uncleared_1772=(|, "cleared_1810", "not_cleared_1810"), uncleared_1810="not_cleared_1810"),
        )),
        # Need a second round with this file as the categories overlap
        atlas_19C_land_use_2 = (filename="Mauritius/Undigitised/atlas_19C_land_use_2.jpg", poly=1, layers=(;
            abandonned=(abandonned_1905_cleared_1968="abdn_1854-1905_cleared_1905-1968",)
        )),
        atlas_19C_land_use = (filename="Mauritius/Undigitised/atlas_19C_land_use.jpg", poly=1, layers=(;
            cleared=(;
                # We assume clearing happened some time before the area became urban
                # and include 1905 urban in 1810 cleared
                cleared_1810=(|, "cleared_1810", "urban_1810", "urban_1905", "cleared_1810_abdn_1905"),
                cleared_1854=(|, "cleared_1810", "urban_1810", "urban_1905", "cleared_1854", "cleared_1810_abdn_1905", "cleared_1854_abdn_1905", "cleared_1854_abdn_1968"),
                cleared_1905=(|, "cleared_1810", "urban_1810", "urban_1905", "cleared_1854", "cleared_1905", "cleared_1905_abdn_1905", "cleared_1854_abdn_1968", "cleared_1905_abdn_1968"),
                cleared_1968=(|, "cleared_1810", "urban_1810", "urban_1905", "cleared_1854", "cleared_1905", "cleared_1968", "cleared_1968_abdn_1968"),
            ),
            lakes = (; lakes="lakes",),
            urban = (; urban_1810=(|, "urban_1810"), urban_1905=(|, "urban_1810", "urban_1905")),
            abandonned=(
                abandonned_1905=(|, "cleared_1810_abdn_1905", "cleared_1854_abdn_1905", "cleared_1905_abdn_1905",),
                abandonned_1968=(|, "cleared_1854_abdn_1968", "cleared_1905_abdn_1905", "cleared_1905_abdn_1968", "cleared_1968_abdn_1968"),
            ),
            # uncleared=(uncleared_1772=(|, "cleared_1810", "not_cleared_1810"), uncleared_1968="not_cleared_1968"),
        )),
        desroches_1773_from_gleadow = (filename="Mauritius/Undigitised/1773_desroches_from_gleadow.jpg", poly=1, layers=(;
            conceded=(conceded_1773="conceded_land"),
        )),
        # fraser_1835_composite_etsy = (filename="Mauritius/Undigitised/1835_fraser_composite_etsy.png", poly=1, layers=(;
        #     cleared=(cleared_1835="cleared",), forest=(forest_1835="uncleared",),
        # )),
        fraser_1835_from_gleadow = (filename="Mauritius/Undigitised/1835_fraser_from_gleadow.jpg", poly=1, layers=(;
            cleared=(cleared_1835="cleared",), forest=(forests_1835="forest",),
        )),
        surveyor_general_1872_from_gleadow = (filename="Mauritius/Undigitised/1872_surveyor_general_from_gleadow.jpg", poly=1, layers=(;
            # cleared=(cleared_1850="cleared_1850-70", cleared_1870=("cleared_1850-70", "cleared_1870-72"), cleared_1872=("cleared_1850-70", "cleared_1870-72")),
            forest=(forest_1850=(|, "forest_1872", "cleared_1850-70", "cleared_1870-72"), forest_1870=(|, "forest_1872", "cleared_1870-72"), forest_1872="forest_1872",),
        )),
     ), reu=(; 
        # cadet_invasives=(filename="Reunion/Undigitised/cadet_invasives.jpg", poly=1, layers=(;)),
        # atlas_vegetation = (filename="Reunion/Undigitised/atlas_vegetation.jpg", poly=1, layers=(;)),
        # atlas_ownership = (filename="Reunion/Undigitised/atlas_ownership.jpg", poly=1, layers=(;)),
        # atlas_1960_population = (filename="Reunion/Undigitised/atlas_1960_population.jpg", poly=1, layers=(;)),
        # "atlas_1960_agriculture" => (filename="Reunion/Undigitised/atlas_agriculture_1960.jpg", poly=1, layers=(;)),
        atlas_1960_agriculture = (filename="Reunion/Undigitised/atlas_agriculture_1960_2.jpg", poly=1, layers=(;
            forest=(|, "rock", "forest", "shrubland", "savannah"),
            agriculture=(|, "cane", "geranium_continuous", "geranium_discontinuous", "tea"),
            reforestation=(|, "casuarina", "acacia", "cryptomeria", "labourdonassia"),
            urban="urban",
        )),
        atlas_1815_agriculture = (filename="Reunion/Undigitised/atlas_1815_agriculture.jpg", poly=1, layers=(;
            uncleared = "forest",
            cleared = (|, "geranium", "vanilla", "cane", "wasteland"),
        )),
        atlas_1780_agriculture = (filename="Reunion/Undigitised/atlas_1780_agriculture.jpg", poly=1, layers=(;
            uncleared = "native",
            # cleared = (!, "forest"),
        )),
        atlas_early_settlement = (filename="Reunion/Undigitised/atlas_early_settlement_cropped.jpg", poly=1, layers=(;
            concessions=(conceded_1715_1765="conceded_1715-1765", conceded_1665_1715="concede_1665-1715"),
        )),
    ))

    files = map(file_details) do island
        map(island) do (; filename, poly, layers)
            (; filename=joinpath(selected_dir, filename), poly, layers)
        end
    end
    return files
end

function make_raster_slices(masks)
    selected_dir = "/home/raf/PhD/Mascarenes/Data/Selected"
    # Copy duplicated file wrap points
    fn = joinpath(selected_dir, "Mauritius/Undigitised/atlas_19C_land_use")
    cp(fn * ".csv", fn * "_2.csv"; force=true)
    files = get_map_files()

    # Load all rasters and make masks
    rasters = map(files) do island_files
        map(island_files) do file
            image_path = file.filename
            raster_path = splitext(image_path)[1] * ".tif"
            raw = view(Raster(raster_path), Band(1))
            json_path = splitext(file.filename)[1] * ".json"
            data = JSON3.read(read(json_path), MapRasterization.MapSelection)
            grouped = map(file.layers) do layer
                _category_raster(raw, data.settings.category_name, layer)
            end
            out = (; raw, grouped)
        end
    end

    # Generate timelines
    mus_timelines = let
        abandonded_1905_cleared_1968 = rasters.mus.atlas_19C_land_use_2.grouped.abandonned.abandonned_1905_cleared_1968
        m = rasters.mus
        cleared_1610 = falses(dims(masks.mus))
        # Atlas dutch period
        cleared_1710 = m.atlas_dutch_period.grouped.cleared .& masks.mus
        # Dutch departure
        cleared_1711 = falses(dims(masks.mus))
        # French arrival
        cleared_1723 = falses(dims(masks.mus))
        # Atlas 18_C + dutch period
        (; cleared_1772, cleared_1810) = map(m.atlas_18C_land_use.grouped.cleared) do A
            A .& masks.mus
        end
        # Atlas 19C + dutch period
        cleared_19C = m.atlas_19C_land_use.grouped.cleared
        (; cleared_1810, cleared_1854, cleared_1905, cleared_1968) = map(m.atlas_19C_land_use.grouped.cleared) do A
            A .& masks.mus
        end
        # Gleadow hand drawn 19C maps
        gleadow_forest = map(m.surveyor_general_1872_from_gleadow.grouped.forest) do A
            A .& masks.mus
        end
        gleadow_cleared_keys = map(keys(gleadow_forest)) do k
            Symbol(replace(string(k), "forest" => "cleared"))
        end
        (; cleared_1850, cleared_1870, cleared_1872) = map(gleadow_forest) do A
            # TODO change this to 1905 after fixing missing areas on 1905 map
            .!(A .| .!(cleared_1968)) .| cleared_1854
        end |> values |> NamedTuple{gleadow_cleared_keys}
        # TODO Delete this after fixing 1905 maps
        cleared_1905 = cleared_1905 .| cleared_1872
        cleared_1992 = m.atlas_1992_agriculture.grouped.cleared .& masks.mus
            # m.atlas_1992_land_use.grouped.cleared

        cleared = (;
            cleared_1610,
            cleared_1710,
            cleared_1711,
            cleared_1723,
            cleared_1772,
            cleared_1810, 
            cleared_1854,
            cleared_1870,
            cleared_1872,
            cleared_1905, 
            cleared_1968,
            cleared_1992,
        )

        abandonned_1610 = falses(dims(masks.mus))
        abandonned_1710 = falses(dims(masks.mus))
        # Dutch departure
        abandonned_1711 = cleared_1710
        # French arrival
        abandonned_1723 = cleared_1710
        abandonned_1772 = abandonned_1723 .& .!(cleared_1772)
        # English arrival
        abandonned_1810 = m.atlas_18C_land_use.grouped.abandonned.abandonned_1810 .& masks.mus
        # 1810 abandonment may be subsequently cleared, so & its inverse
        abandonned_1854 = abandonned_1810 .& .!(cleared_1854)
        abandonned_1870 = abandonned_1810 .& .!(cleared_1870)
        abandonned_1872 = abandonned_1810 .& .!(cleared_1872)
        abandonned_1905 = m.atlas_19C_land_use_2.grouped.abandonned.abandonned_1905_cleared_1968 .|
            m.atlas_19C_land_use.grouped.abandonned.abandonned_1905 .|
            (abandonned_1810 .& .!(cleared.cleared_1905)) .&
            masks.mus
        abandonned_1968 = abandonned_1810 .|
            m.atlas_19C_land_use.grouped.abandonned.abandonned_1968 .|
            m.atlas_19C_land_use.grouped.abandonned.abandonned_1905 .|
            (abandonned_1810 .& .!(cleared_1968)) .&
            masks.mus
        abandonned_1992 = ((abandonned_1810 .| abandonned_1968) .|
            # There is no abandonment data for 1992 so we use the difference with 1968
            (cleared_1968 .& .!(cleared_1992))) .& 
            # And remove anything cleared in 1992
            .!(cleared_1992)

        abandonned = (;
            abandonned_1610,  
            abandonned_1710,
            abandonned_1711,
            abandonned_1723,
            abandonned_1772,
            abandonned_1810,
            abandonned_1854,
            abandonned_1870,
            abandonned_1872,
            abandonned_1905,
            abandonned_1968,
            abandonned_1992,
        )

        urban = (;
            urban_1763=m.atlas_18C_land_use.grouped.urban.urban_1763,
            # urban_1810=m.atlas_18C_land_use.grouped.urban.urban_1810,
            urban_1810=m.atlas_19C_land_use.grouped.urban.urban_1810,
            urban_1905=m.atlas_19C_land_use.grouped.urban.urban_1905,
            urban_1992=m.atlas_1992_agriculture.grouped.urban,
        )
        # Multiply by 2 so that abandonned can be 1
        cleared2 = map(A -> 2 .* A, cleared)
        # Add abandonned values where they exist
        _combine(A, B) = broadcast((a, b) -> b ? 1 : a, A, B)
        # Combine cleared and abandonned Bool masks into Int rasters
        combined_keys = map(x -> Symbol("combined" * string(x)[end-4:end]), keys(cleared))
        combined = NamedTuple{combined_keys}(map(_combine, values(cleared2), values(abandonned)))

        (; cleared, abandonned, urban, combined)
    end

    reu_timelines = let
        r = rasters.reu
        cleared_1780 = (.!(r.atlas_1780_agriculture.grouped.uncleared)) .& masks.reu
        cleared_1715 = r.atlas_early_settlement.grouped.concessions.conceded_1665_1715 .& cleared_1780 .& masks.reu
        cleared_1765 = ((r.atlas_early_settlement.grouped.concessions.conceded_1715_1765 .| cleared_1715) .& cleared_1780) .& masks.reu
        cleared_1815 = (r.atlas_1815_agriculture.grouped.cleared .| cleared_1780) .& masks.reu
        cleared_1960 = (.!(r.atlas_1960_agriculture.grouped.forest) .| cleared_1815) .& masks.reu
        cleared = (; cleared_1715, cleared_1765, cleared_1780, cleared_1815, cleared_1960,)
        # We don't have abandonment data for Reunion
        (; cleared)
    end

    # output NamedTuples with (:raw, :grouped, :timeline) keys
    return (;
        mus=(files=rasters.mus, timelines=mus_timelines),
        reu=(files=rasters.reu, timelines=reu_timelines),
    )
end # make_raster_slices

function _category_raster(raster::Raster, category_names::Vector, layers::NamedTuple)
    map(layers) do layer
        _category_raster(raster, category_names, layer)
    end
end
function _category_raster(raster::Raster, category_names::Vector, layer_components::Tuple)
    # The first argument is a function, run it on the rest of the slices
    op = first(layer_components)
    layers = map(l -> _category_raster(raster, category_names, l), Base.tail(layer_components))
    broadcast(op, layers...)
end
function _category_raster(raster, category_names::Vector, category::String)
    I = findall(==(category), map(String, category_names))
    if length(I) == 0
        error("could not find $category in $(category_names)")
    end
    out = raster .== first(I)
    foreach(I[2:end]) do i
        out .|= raster .== first(i)
    end
    return out
end
function _category_to_raster(raster::Raster, category_names::Vector, x)
    error("slice must be a String, NamedTuple or Tuple, got $x")
end


open_output(T, x::NamedTuple) = open_output(x.filename)
function open_output(T, filename::String)
    json_path = splitext(filename)[1] * ".json"
    return isfile(json_path) ? JSON3.read(read(json_path), T) : nothing
end

open_warp_points(x::NamedTuple) = open_warp_points(x.filename)
function open_warp_points(filename::String)
    csv_path = splitext(filename)[1] * ".csv"
    return isfile(csv_path) ? CSV.read(csv_path, DataFrame) : nothing
end

function warp_to_raster(img_path::String, template::Raster; object_type=MapRasterization.MapSelection, edit=false, kw...)
    img = load_image(img_path)
    csv_path = splitext(img_path)[1] * ".csv"
    points = isfile(csv_path) ? open_warp_points(img_path) : nothing
    if edit || !isfile(csv_path)
        warp_points = if isnothing(points)
            MapRasterization.select_warp_points(img;
                template=reverse(template; dims=Y()), kw...
            )
        else
            :x_known in names(points) && rename!(points, [:x_known => :x_a, :y_known => :y_a, :x_unknown => :x_b, :y_unknown => :y_b])
            # if :x_a in names(points)
                MapRasterization.select_warp_points(img;
                    template=reverse(template; dims=Y()), points, kw...
                )
            # else
                # MapRasterization.select_warp_points(Float64.(Gray.(img));
                    # template=reverse(template; dims=Y()), missingval=missingval(template),
                # )
            # end
        end
        df = DataFrame(warp_points)
        CSV.write(csv_path, df)
    end
    if isfile(splitext(img_path)[1] * ".json")
        df = CSV.read(csv_path, DataFrame)
        output = open_output(object_type, img_path)
        if object_type <: MapRasterization.MapSelection
            out = Int.(reshape(output.output, size(img)))
            rs = MapRasterization.applywarp(out; template, points=df, missingval=0)
            raster_path = splitext(img_path)[1] * ".tif"
            write(raster_path, rs)
            rs = mask(Raster(raster_path); with=template)
            write(raster_path, rs)
            display(Plots.plot(rs))
            return rs
        elseif object_type <: MapRasterization.CategoryShapes
            warped_geoms = map(output.shapes) do sh
                geoms = Polygon.(sh)
                warped = MapRasterization.applywarp(geoms; template, points=df)
                map(warped) do g
                    collect(GeoInterface.getpoint(g))
                end
            end
            warped = MapRasterization.CategoryShapes{Polygon}(warped_geoms, output.category_names)
            # GeoJSON.write(splitext(img_path)[1] * "_warped.geojson", warped)
            return warped
        end
    end
end

function choose_categories(img_path::String; output=open_output(MapRasterization.MapSelection, img_path), save=true)
    img = load_image(img_path)
    if isnothing(output)
        output = MapRasterization.selectcolors(img; ncategories=5)
    else
        output = MapRasterization.selectcolors(img, output)
    end
    if save
        json_path = splitext(img_path)[1] * ".json"
        if isfile(json_path)
            backup_path = splitext(img_path)[1] * "_backup.json"
            cp(json_path, backup_path; force=true)
        end
        write(json_path, JSON3.write(output))
    end
    return output
end

load_image(img_path::String) = RGB{Float64}.(load(img_path) |> rotr90)
