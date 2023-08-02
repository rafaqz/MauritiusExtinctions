using JSON3, MapRasterization, GeoInterface, Rasters, FileIO, ImageIO, DataFrames, CSV, GeoJSON, Colors
using DimensionalData.LookupArrays


function define_map_files(; path = "/home/raf/PhD/Mascarenes/Data/Selected")
    # Here we define:
    # - all of the files we use
    # - what the land-cover categories are called for the file
    # - a map between those categories in our main land cover categories
    #
    # These use either a single String or Tuples of strings staring with a function.
    # The function will later be broadcasted over masks of the separate layers to combine them
    # Mostly is `|` which is "or" so we make a mask of values that are true in one or the other file
    file_details = (mus=(;
        atlas_18C_land_use = (filename="Mauritius/Undigitised/atlas_18C_land_use.jpg", poly=1, layers=(;
            cleared=(
                cleared_1772=["urban_1763", "cleared_1772", "abandoned_1810"],
                cleared_1810=["urban_1763", "cleared_1772", "urban_1810", "cleared_1810"],
            ),
            urban = (urban_1763="urban_1763", urban_1810=["urban_1763", "urban_1810"]),
            abandoned=(abandoned_1810="abandoned_1810",),
            # uncleared=(uncleared_1772=["cleared_1810", "not_cleared_1810"), uncleared_1810="not_cleared_1810"),
        )),
        atlas_1992_vegetation = (filename="Mauritius/Undigitised/atlas_1992_vegetation.jpg", poly=1, layers=(;)),
        atlas_1992_agriculture = (filename="Mauritius/Undigitised/atlas_1992_agriculture.jpg", poly=1, layers=(;
            urban="urban",
            cleared=["cane", "forage", "tea", "market_gardens", "cane_or_market_gardens", 
                     "cane_or_tea", "tea_or_market_gardens", "cane_or_fruit", "pasture"],
            forestry="forestry",
            water="lakes",
            unsure="forest",
        )),
        atlas_1992_land_use = (filename="Mauritius/Undigitised/atlas_1992_land_use.jpg", poly=1, layers=(;
            urban=["urban", "other_state_land_urban"],
            cleared=["small_properties", "medium_properties", "large_properties", "rose_bell"],
            unsure= ["other_state_land", "forest", "mountain_reserves", "tea_development_forest", 
                     "private_forest_or_wasteland"],
        )),
        vegetation = (filename="Mauritius/Undigitised/page33_mauritius_vegetation_colored.png", poly=1, layers=(;)),
        atlas_dutch_period = (filename="Mauritius/Undigitised/atlas_dutch_period.jpg", poly=1, layers=(;
            uncleared="undisturbed",
            ebony_harvest="ebony_harvest",
            cleared="cleared",
        )),
        atlas_19C_land_use = (filename="Mauritius/Undigitised/atlas_19C_land_use.jpg", poly=1, layers=(;
            cleared=(;
                # We assume that clearing happened some time before the area
                # became urban, so we include 1905 urban in 1810 cleared
                # because these urban areas of the map hide information about clearing.
                cleared_1810=["cleared_1810", "urban_1810", "urban_1905", "cleared_1810_abdn_1905"],
                cleared_1854=["cleared_1810", "urban_1810", "urban_1905", "cleared_1854", "cleared_1810_abdn_1905", "cleared_1854_abdn_1905", "cleared_1854_abdn_1968"],
                cleared_1905=["cleared_1810", "urban_1810", "urban_1905", "cleared_1854", "cleared_1905", "cleared_1854_abdn_1968", "cleared_1905_abdn_1968"],
                cleared_1968=["cleared_1810", "urban_1810", "urban_1905", "cleared_1854", "cleared_1905", "cleared_1968"],
            ),
            water = "lakes",
            urban = (; urban_1810="urban_1810", urban_1905=["urban_1810", "urban_1905"]),
            abandoned=(
                abandoned_1905=["cleared_1810_abdn_1905", "cleared_1854_abdn_1905", "cleared_1905_abdn_1905",],
                abandoned_1968=["cleared_1810_abdn_1905", "cleared_1854_abdn_1968", "cleared_1905_abdn_1905", "cleared_1905_abdn_1968", "cleared_1968_abdn_1968"],
            ),
            # uncleared=(uncleared_1772=["cleared_1810", "not_cleared_1810"), uncleared_1968="not_cleared_1968"),
        )),
        # We need a second round with this file as the categories overlap 
        atlas_19C_land_use_2 = (filename="Mauritius/Undigitised/atlas_19C_land_use_2.jpg", poly=1, layers=(;
            abandoned=(abandoned_1905_cleared_1968="abdn_1854-1905_cleared_1905-1968",),
            # urban = (; urban_1968=["urban_1968")),
        )),
        desroches_1773_from_gleadow = (filename="Mauritius/Undigitised/1773_desroches_from_gleadow.jpg", poly=1, layers=(;
            conceded=(conceded_1773="conceded_land"),
        )),
        # fraser_1835_composite_etsy = (filename="Mauritius/Undigitised/1835_fraser_composite_etsy.png", poly=1, layers=(;
        #     cleared=(cleared_1835="cleared",), uncleared=(forest_1835="uncleared",),
        # )),
        fraser_1835_from_gleadow = (filename="Mauritius/Undigitised/1835_fraser_from_gleadow.jpg", poly=1, layers=(;
            cleared=(cleared_1835="sea",), # sea and cleared are swapped
            uncleared=(forest_1835="forest",),
        )),
        surveyor_general_1872_from_gleadow = (filename="Mauritius/Undigitised/1872_surveyor_general_from_gleadow.jpg", poly=1, layers=(;
            # cleared=(cleared_1850="cleared_1850-70", cleared_1870=("cleared_1850-70", "cleared_1870-72"), cleared_1872=("cleared_1850-70", "cleared_1870-72")),
            forest=(forest_1850=["forest_1872", "cleared_1850-70", "cleared_1870-72"],
                    forest_1870=["forest_1872", "cleared_1870-72"],
                    forest_1872="forest_1872",
            ),
        )),
        landcover_1965 = (filename="Mauritius/Undigitised/mus_landuse_1965_100_hi_c.pdf", poly=1, layers=(;
            urban="Built_up",
            cleared=["Sugar", "Vegetables", "Tea"],
            forestry="Forest_plantation",
            water="Reservoirs",
            unsure=["Forest_natural", "Rock", "Swamps", "Scrub", "Savannah"], 
        )),
     ), reu=(;
        # cadet_invasives=(filename="Reunion/Undigitised/cadet_invasives.jpg", poly=1, layers=(;)),
        # atlas_vegetation = (filename="Reunion/Undigitised/atlas_vegetation.jpg", poly=1, layers=(;)),
        # atlas_ownership = (filename="Reunion/Undigitised/atlas_ownership.jpg", poly=1, layers=(;)),
        # atlas_1960_population = (filename="Reunion/Undigitised/atlas_1960_population.jpg", poly=1, layers=(;)),
        # "atlas_1960_agriculture" => (filename="Reunion/Undigitised/atlas_agriculture_1960.jpg", poly=1, layers=(;)),
        # atlas_population_1967 = (filename="Reunion/Undigitised/atlas_population_1967.jpg", poly=1, layers=(;)),
        atlas_1960_agriculture = (filename="Reunion/Undigitised/atlas_agriculture_1960_2.jpg", poly=1, layers=(;
            uncleared=["rock", "forest", "shrubland", "savannah"],
            agriculture=["cane", "geranium_continuous", "geranium_discontinuous", "tea"],
            forestry=["casuarina", "acacia", "cryptomeria", "labourdonassia"],
            urban="urban",
        )),
        atlas_1815_agriculture = (filename="Reunion/Undigitised/atlas_1815_agriculture.jpg", poly=1, layers=(;
            uncleared = "forest",
            cleared = ["geranium", "vanilla", "cane"],
            abandonned = "wasteland"
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
            (; filename=joinpath(path, filename), poly, layers)
        end
    end
    return files
end

function make_raster_slices(masks, categories; path="/home/raf/PhD/Mascarenes/Data/Selected")
    # Copy duplicated file wrap points
    fn = joinpath(path, "Mauritius/Undigitised/atlas_19C_land_use")
    cp(fn * ".csv", fn * "_2.csv"; force=true)
    files = define_map_files()

    # Load all rasters and make masks layers for all land-cover classes.
    # We get the numbers form the saved ".json" file for the project.
    rasters = map(files) do island_files
        map(island_files) do file
            image_path = file.filename
            raster_path = splitext(image_path)[1] * ".tif"
            raw = Raster(raster_path)
            json_path = splitext(file.filename)[1] * ".json"
            data = JSON3.read(read(json_path), MapRasterization.MapSelection)
            grouped = map(file.layers) do layer
                _category_raster(raw, data.settings.category_name, layer)
            end
            out = (; raw, grouped)
        end
    end

    # Generate timelines
    # Here we combined slices from multiple separate maps to 
    # make time slices for each land-cover type.
    mus_timelines = let
        m = rasters.mus
        emptymask = falses(dims(masks.mus))

        forestry_1965 = m.landcover_1965.grouped.forestry
        forestry_1992 = m.atlas_1992_agriculture.grouped.forestry .| forestry_1965

        # Assume we never abandon urban areas
        urban_1600 = emptymask
        urban_1721 = emptymask
        urban_1763 = m.atlas_18C_land_use.grouped.urban.urban_1763
        urban_1810 = m.atlas_19C_land_use.grouped.urban.urban_1810 .| urban_1763
        urban_1905 = m.atlas_19C_land_use.grouped.urban.urban_1905 .| urban_1810
        urban_1965 = m.landcover_1965.grouped.urban .| urban_1905
        urban_1992 = m.atlas_1992_agriculture.grouped.urban .| urban_1965

        abandonded_1905_cleared_1968 = m.atlas_19C_land_use_2.grouped.abandoned.abandoned_1905_cleared_1968

        # Clearing
        # No clearing has occurred yet
        cleared_1600 = emptymask
        cleared_1638 = emptymask
        # Atlas dutch period
        cleared_1709 = m.atlas_dutch_period.grouped.cleared
        # Dutch departure abandonment
        cleared_1711 = emptymask
        # French arrival still abandoned
        cleared_1721 = emptymask
        # Atlas 18_C
        cleared_1772 = m.atlas_18C_land_use.grouped.cleared.cleared_1772
        cleared_1810 = m.atlas_18C_land_use.grouped.cleared.cleared_1810
        # Atlas 19C
        cleared_19C = m.atlas_19C_land_use.grouped.cleared
        cleared_1810 = cleared_19C.cleared_1810
        cleared_1854 = cleared_19C.cleared_1854
        cleared_1905 = cleared_19C.cleared_1905
        cleared_1965 = m.landcover_1965.grouped.cleared
        cleared_1968 = cleared_19C.cleared_1968 .| cleared_1965

        fraser_forest = m.fraser_1835_from_gleadow.grouped.uncleared.forest_1835
        # Gleadow hand drawn 19C maps
        # We try to reduce cartographical innacuracies of the projection by masking
        # the gleadow and fraser maps with the previous and subsequent atlas maps.
        # This will have the effect of reducing the total cleared area for thes
        # frames by a small amount
        # TODO Delete this after fixing 1905 maps?
        # cleared_1905 = cleared_1905 .| cleared_1872
        cleared_1992 = m.atlas_1992_agriculture.grouped.cleared
            # m.atlas_1992_land_use.grouped.cleared

        # Abandoned cleared land
        abandoned_1600 = emptymask
        abandoned_1638 = emptymask
        abandoned_1709 = emptymask
        # Dutch departure
        abandoned_1711 = cleared_1709
        # French arrival
        abandoned_1721 = abandoned_1711
        abandoned_1772 = (abandoned_1721 .| cleared_1721 .| urban_1721) .& .!(cleared_1772 .| urban_1763)
        # English arrival
        abandoned_1810 = (((abandoned_1772 .| cleared_1772 .| urban_1810) .& .!(cleared_1810 .| urban_1810)) .|
                          m.atlas_18C_land_use.grouped.abandoned.abandoned_1810)
        # 1810 abandonment may be subsequently cleared, so & its inverse
        # We don't have urban data for this time, but the cleared maps appear to include urban
        abandoned_1854 = ((abandoned_1810 .| cleared_1810) .& .!(cleared_1854 .| urban_1810))
        # abandoned_1870 = ((abandoned_1854 .| cleared_1854) .& .!(cleared_1870 .| urban_1810))
        # abandoned_1872 = ((abandoned_1870 .| cleared_1870) .& .!(cleared_1872 .| urban_1810))
        abandoned_1905 = m.atlas_19C_land_use_2.grouped.abandoned.abandoned_1905_cleared_1968 .|
            m.atlas_19C_land_use.grouped.abandoned.abandoned_1905 .& .!(cleared_1905)
        abandoned_1965 = ((abandoned_1905 .| cleared_1905) .& .!(cleared_1965))
        abandoned_1968 = m.atlas_19C_land_use.grouped.abandoned.abandoned_1968 .| abandoned_1965
        # There is no abandonment data for 1992 so we use 
        # previous abandonned land and the difference with 1968 cleared land
        abandoned_1992 = (abandoned_1810 .| abandoned_1968 .| cleared_1968) .& 
            # And remove anything cleared or urban in 1992
            .!(cleared_1992 .| urban_1992 .| forestry_1992)

        # We need abandonement data from 1810 so this is moved later
        # Cleared areas in 1835 must not be uncleared again later in 1854
        cleared_1835 = .!(fraser_forest) .& (cleared_1854 .| abandoned_1854)
        # Uncleared areas that were previously cleared or abandoned
        # in 1810 or abandoned in 1854 are also abandoned in 1835
        abandoned_1835 = fraser_forest .& (abandoned_1810 .| cleared_1810 .| abandoned_1854)

        sg_forest = m.surveyor_general_1872_from_gleadow.grouped.forest
        uncleared_1905 = .!(urban_1905 .| cleared_1905 .| abandoned_1905)
        cleared_1850 = .!(uncleared_1905 .& sg_forest.forest_1850) .& .!(urban_1810)
        cleared_1870 = .!(uncleared_1905 .& sg_forest.forest_1870) .| cleared_1854 .& .!(urban_1810)
        cleared_1872 = .!(uncleared_1905 .& sg_forest.forest_1872) .| cleared_1854 .& .!(urban_1810)

        cleared = [
            1600=>masks.mus .& cleared_1600,
            1638=>masks.mus .& cleared_1638,
            1709=>masks.mus .& cleared_1709,
            1711=>masks.mus .& cleared_1711,
            1721=>masks.mus .& cleared_1721,
            1772=>masks.mus .& cleared_1772,
            1810=>masks.mus .& cleared_1810,
            1835=>masks.mus .& cleared_1835,
            # 1850=>masks.mus .& cleared_1850,
            1854=>masks.mus .& cleared_1854,
            # 1870=>masks.mus .& cleared_1870,
            # 1872=>masks.mus .& cleared_1872,
            1905=>masks.mus .& cleared_1905,
            1965=>masks.mus .& cleared_1965,
            1968=>masks.mus .& cleared_1968,
            1992=>masks.mus .& cleared_1992,
        ]

        cleared = RasterSeries(last.(cleared), Ti(first.(cleared); 
            span=Irregular(1500, first(last(cleared))),
            sampling=Intervals(End()),
        ))
        abandoned = [
            1600=>masks.mus .& abandoned_1600,
            1638=>masks.mus .& abandoned_1600,
            1709=>masks.mus .& abandoned_1709,
            1711=>masks.mus .& abandoned_1711,
            1721=>masks.mus .& abandoned_1721,
            1772=>masks.mus .& abandoned_1772,
            1810=>masks.mus .& abandoned_1810,
            1835=>masks.mus .& abandoned_1835,
            1854=>abandoned_1854,
            # 1870=>abandoned_1870,
            # 1872=>abandoned_1872,
            1905=>masks.mus .& abandoned_1905,
            1968=>masks.mus .& abandoned_1968,
            1992=>masks.mus .& abandoned_1992,
        ]
        abandoned = RasterSeries(last.(abandoned), Ti(first.(abandoned); 
            span=Irregular(1500, first(last(abandoned))),
            sampling=Intervals(End())
        ))
        urban = [
            1600=>masks.mus .& urban_1600,
            1721=>masks.mus .& urban_1721,
            1763=>masks.mus .& urban_1763,
            1810=>masks.mus .& urban_1810,
            1905=>masks.mus .& urban_1905,
            1965=>masks.mus .& urban_1965,
            1992=>masks.mus .& urban_1992,
        ]
        urban = RasterSeries(last.(urban), Ti(first.(urban); 
            span=Irregular(1500, first(last(urban))),
            sampling=Intervals(End()),
        ))
        forestry = [
            1600=>emptymask,
            1905=>emptymask,
            1965=>masks.mus .& forestry_1965,
            1992=>masks.mus .& forestry_1992,
        ]
        forestry = RasterSeries(last.(forestry), Ti(first.(forestry); 
            span=Irregular(1500, first(last(forestry))),
            sampling=Intervals(End()),
        ))

        water_1905 = m.atlas_19C_land_use.grouped.water
        water_1965 = m.landcover_1965.grouped.water .| water_1905
        water_1992 = m.atlas_1992_agriculture.grouped.water .| water_1965
        water = [
            1905=>masks.mus .& water_1905
            1965=>masks.mus .& water_1965
            1992=>masks.mus .& water_1992
        ]
        water = RasterSeries(last.(water), Ti(first.(water); 
            span=Irregular(1500, first(last(water))),
            sampling=Intervals(End()),
        ))
        function _combine(W, F, C, A, U)
            A = broadcast(W, F, C, A, U, masks.mus) do w, f, c, a, u, m
                if !m
                    0
                elseif w
                    categories.water
                elseif f
                    categories.forestry
                elseif u
                    categories.urban
                elseif a
                    categories.abandoned
                elseif c
                    categories.cleared
                else
                    categories.native
                end
            end
            rebuild(A; missingval=0)
        end
        # Combine cleared and abandoned Bool masks into Int rasters
        lc = RasterSeries(map(lookup(cleared, Ti)) do t
            i = Ti(Contains(t))
            _combine(water[i], forestry[i], cleared[i], abandoned[i], urban[i])
        end, dims(cleared, Ti))
        (; cleared, abandoned, urban, forestry, water, lc)
    end

    reu_timelines = let
        r = rasters.reu
        cleared_1780 = (.!(r.atlas_1780_agriculture.grouped.uncleared)) .& masks.reu
        cleared_1715 = r.atlas_early_settlement.grouped.concessions.conceded_1665_1715 .& cleared_1780 .& masks.reu
        cleared_1765 = ((r.atlas_early_settlement.grouped.concessions.conceded_1715_1765 .| cleared_1715) .& cleared_1780) .& masks.reu
        cleared_1815 = (r.atlas_1815_agriculture.grouped.cleared .| cleared_1780) .& masks.reu
        cleared_1960 = (.!(r.atlas_1960_agriculture.grouped.uncleared) .| cleared_1815) .& masks.reu
        # cleared = Rasters.combine(RasterSeries(cleared, Ti([1715, 1765, 1780, 1815, 1960]; sampling=Intervals(End()))))
        cleared = RasterSeries(
            [cleared_1715, cleared_1765, cleared_1780, cleared_1815, cleared_1960],
            Ti([1715, 1765, 1780, 1815, 1960]; sampling=Intervals(End())),
        )
        urban = RasterSeries(
            [cleared_1715, cleared_1765, cleared_1780, cleared_1815, cleared_1960],
            Ti([1715, 1765, 1780, 1815, 1960]; sampling=Intervals(End())),
        )
        # We don't have abandonment data for Reunion
        (; urban, cleared)
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
function _category_raster(raster::Raster, category_names::Vector, layer_components::Vector{String})
    layers = map(l -> _category_raster(raster, category_names, l), layer_components)
    broadcast(|, layers...)
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
    error("slice must be a String, NamedTuple or Vector{String}, got $x")
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

function warp_to_raster(img_path::String, template::Raster;
    object_type=MapRasterization.MapSelection, edit=false, save=true, kw...
)
    img = load_image(img_path)
    csv_path = splitext(img_path)[1] * ".csv"
    points = isfile(csv_path) ? open_warp_points(img_path) : nothing
    if edit || !isfile(csv_path)
        warp_points = if isnothing(points)
            MapRasterization.click_warp(img;
                template=reverse(template; dims=Y()), kw...
            )
        else
            :x_known in names(points) && rename!(points, [:x_known => :x_a, :y_known => :y_a, :x_unknown => :x_b, :y_unknown => :y_b])
            # if :x_a in names(points)
                MapRasterization.click_warp(img;
                    template=reverse(template; dims=Y()), points, kw...
                )
            # else
                # MapRasterization.click_warp(Float64.(Gray.(img));
                    # template=reverse(template; dims=Y()), missingval=missingval(template),
                # )
            # end
        end
        if save
            df = DataFrame(warp_points)
            CSV.write(csv_path, df)
        end
    end
    if save && isfile(splitext(img_path)[1] * ".json")
        df = CSV.read(csv_path, DataFrame)
        output = open_output(object_type, img_path)
        poly = 1
        if object_type <: MapRasterization.MapSelection
            out = Int.(reshape(output.output, size(img)))
            warper = MapRasterization.Warper(df, template, poly)
            rs = MapRasterization.warp(warper, out; missingval=0)
            raster_path = splitext(img_path)[1] * ".tif"
            # write(raster_path, rs)
            rs = mask(Raster(raster_path); with=template)
            write(raster_path, rs)
            return rs
        elseif object_type <: MapRasterization.CategoryShapes
            warper = MapRasterization.Warper(df, template, poly)
            warped_geoms = map(output.shapes) do sh
                geoms = Polygon.(sh)
                warped = MapRasterization.warp(warper, geoms)
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

function choose_categories(img_path::String;
    save=true, restart=false,
    output=restart ? nothing : open_output(MapRasterization.MapSelection, img_path),
)
    img = load_image(img_path)
    if isnothing(output)
        cs = MapRasterization.CategorySelector(img)
    else
        cs = MapRasterization.CategorySelector(img, output)
    end
    if save
        output = MapRasterization.MapSelection(cs)
        json_path = splitext(img_path)[1] * ".json"
        if isfile(json_path)
            backup_path = splitext(img_path)[1] * "_backup.json"
            cp(json_path, backup_path; force=true)
        end
        write(json_path, JSON3.write(output))
    end
    return cs
end

load_image(img_path::String) = RGB{Float64}.(load(img_path) |> rotr90)
