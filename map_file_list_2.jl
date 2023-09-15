using JSON3, MapRasterization, GeoInterface, Rasters, FileIO, ImageIO, DataFrames, CSV, GeoJSON, Colors
using DimensionalData.LookupArrays
using InvertedIndices, Setfield

function define_map_files(; path = "/home/raf/PhD/Mascarenes")
    @show path
    # Here we define:
    # - all of the files we use
    # - what the land-cover categories are called for the file
    # - a map between those categories in our main land cover categories
    #
    # These use either a single String or Tuples of strings staring with a function.
    # The function will later be broadcasted over masks of the separate layers to combine them
    # Mostly is `|` which is "or" so we make a mask of values that are true in one or the other file
    file_details = (mus=(;
        atlas_dutch_period = "$path/Data/Selected/Mauritius/Undigitised/atlas_dutch_period.jpg" => (;
            native=[
                1600 => :mask,
                1709 => ["ebony_harvest", "undisturbed"],
                1710 => ["ebony_harvest", "undisturbed"],
            ],
            # undisturbed=1709 => "undisturbed",
            # disturbed=1709 => "ebony_harvest",
            cleared=1709 => "cleared",
            abandoned=1710 => "cleared",
        ),
        atlas_18C_land_use = "$path/Data/Selected/Mauritius/Undigitised/atlas_18C_land_use.jpg" => (;
            native=[
                1763 => ["not_cleared_1810", "cleared_1772", "cleared_1810"],
                1772 => ["not_cleared_1810", "cleared_1810"],
                1810 => ["not_cleared_1810"],
            ],
            cleared=[
                1763 => ["cleared_1772", "abandoned_1810", "urban_1810"],
                1772 => ["cleared_1772", "abandoned_1810", "urban_1810"],
                1810 => ["cleared_1772", "cleared_1810"],
            ],
            urban = [
                1763 => "urban_1763",
                1772 => ["urban_1763", "urban_1810"],
                1810 => ["urban_1763", "urban_1810"],
            ],
            abandoned=[
                1763 => "abandoned_1810",
                1772 => "abandoned_1810",
                1810 => "abandoned_1810",
            ],
        ),
        fraser_1835_from_gleadow = "$path/Data/Selected/Mauritius/Undigitised/1835_fraser_from_gleadow.jpg" => (;
            cleared=1835 => "sea", # sea and cleared are swapped
            urban=1835 => "sea", # we don't know what part of the cleared area was urban
            abandoned=1835 => "sea", # or abandoned
            native=1835 => "forest",
            # water=1835 => "forest",
        ),
        atlas_19C_land_use = "$path/Data/Selected/Mauritius/Undigitised/atlas_19C_land_use.jpg" => (;
            native=[
                1810 => ["not_cleared_1968", "cleared_1854", "cleared_1854_abdn_1905", "cleared_1854_abdn_1968", "cleared_1905", "cleared_1968", "cleared_1968_abdn_1968"],
                1854 => ["not_cleared_1968", "cleared_1905", "cleared_1968", "cleared_1968_abdn_1968"],
                1905 => ["not_cleared_1968", "cleared_1968", "cleared_1968_abdn_1968"], 
                # 1968 => "not_cleared_1968",
            ],
            cleared=[
                # We assume that clearing happened some time before the area
                # became urban, so we include 1905 urban in 1810 cleared
                # because these urban areas of the map hide information about clearing.
                # 1810 => ["cleared_1810", "urban_1905", "cleared_1810_abdn_1905"],
                1854 => ["cleared_1810", "cleared_1854", "urban_1905", "cleared_1810_abdn_1905", "cleared_1810_abdn_1968", "cleared_1854_abdn_1905", "cleared_1854_abdn_1968"],
                1905 => ["cleared_1810", "cleared_1854", "cleared_1905", "cleared_1810_abdn_1968", "cleared_1854_abdn_1968", "cleared_1905_abdn_1968"],
                # This is worse than the 1965 map
                # 1968 => ["cleared_1810", "cleared_1854", "cleared_1905", "cleared_1968"],

            ],
            abandoned=[
                1854 => ["cleared_1810_abdn_1905"],
                1905 => ["cleared_1810_abdn_1905", "cleared_1854_abdn_1905", "cleared_1905_abdn_1905",],
                # 1968 => ["cleared_1810_abdn_1905", "cleared_1854_abdn_1905", "cleared_1905_abdn_1905", "cleared_1810_abdn_1968", "cleared_1854_abdn_1968", "cleared_1905_abdn_1968", "cleared_1968_abdn_1968"],
            ],
            urban=[
                1810 => "urban_1810",
                1854 => ["urban_1810", "urban_1905"],
                1905 => ["urban_1810", "urban_1905"],
                # 1968 => ["urban_1810", "urban_1905", "cleared_1810", "cleared_1854", "cleared_1905"],
            ],
            forestry = [
                # 1968 => ["not_cleared_1968", "cleared_1810_abdn_1905", "cleared_1854_abdn_1905", "cleared_1905_abdn_1905", "cleared_1810_abdn_1968", "cleared_1854_abdn_1968", "cleared_1905_abdn_1968", "cleared_1968_abdn_1968"],
            ],
            water=[
                1810 => "lakes",
                1854 => "lakes",
                1905 => "lakes",
                # 1968 => "lakes",
            ],
        ),
        landcover_1965 = "$path/Data/Selected/Mauritius/Undigitised/mus_landuse_1965_100_hi_c.pdf" => (;
            native=1965 => ["Forest_natural", "Swamps", "Rock", "Scrub", "Savannah"],
            cleared=1965 => ["Sugar", "Vegetables", "Tea"],
            abandoned=1965 => ["Rock", "Scrub", "Savannah"],
            urban=1965 => "Built_up",
            forestry=1965 => "Forest_plantation",
            water=[
                1965 => "Reservoirs",
                2020 => "Reservoirs",
            ],
        ),
        # atlas_1992_vegetation = "$path/Data/Selected/Mauritius/Undigitised/atlas_1992_vegetation.jpg" => (;),
        atlas_1992_agriculture = "$path/Data/Selected/Mauritius/Undigitised/atlas_1992_agriculture.jpg" =>
        (;
            native=1992 => "forest",
            cleared=1992 => [
                "cane", "forage", "tea", "market_gardens", "cane_or_market_gardens",
                "cane_or_tea", "tea_or_market_gardens", "cane_or_fruit", "pasture",
            ],
            urban=1992 => "urban",
            forestry=[
                1992 => "forestry",
                2020 => "forestry",
                2021 => "forestry",
            ],
            water=1992 => "lakes",
            abandoned=1992 => ["pasture", "forest"],
        ),
        # We need a second round with this file as the categories overlap
        # This will overwrite anything incorrect in the first file for the specified year
        atlas_19C_land_use_2 = "$path/Data/Selected/Mauritius/Undigitised/atlas_19C_land_use_2.jpg" => (;
            abandoned=[
                1905 => "abdn_1854-1905_cleared_1905-1968",
            ],
            cleared=[
                1854 => "abdn_1854-1905_cleared_1905-1968",
                # 1968 => "abdn_1854-1905_cleared_1905-1968",
            ],
        ),
        # desroches_1773_from_gleadow = "$path/Data/Selected/Mauritius/Undigitised/1773_desroches_from_gleadow.jpg" => (;
        #     conceded=1773 => "conceded_land",
        # ),
        # atlas_1992_land_use = "$path/Data/Selected/Mauritius/Undigitised/atlas_1992_land_use.jpg" => (;
        #     urban=1992 => ["urban", "other_state_land_urban"],
        #     cleared=1992 => [
        #         "small_properties", "medium_properties", "large_properties", "rose_bell",
        #         "other_state_land", "mountain_reserves", "tea_development_forest",
        #         "forest", "private_forest_or_wasteland",
        #     ],
        #     abandoned=1992 => [
        #         "small_properties", "medium_properties", "large_properties", "rose_bell",
        #         "other_state_land", "mountain_reserves", "tea_development_forest",
        #         "forest", "private_forest_or_wasteland",
        #     ],
        #     forestry=[
        #         1992 => [
        #             "other_state_land", "mountain_reserves", "tea_development_forest",
        #             "forest", "private_forest_or_wasteland",
        #         ],
        #         # Force forestry to stop growing after 1992 we know it stopped.
        #         # 2021 => [
        #         #     "other_state_land", "mountain_reserves", "tea_development_forest",
        #         #     "forest", "private_forest_or_wasteland",
        #         # ],
        #     ],
        #     native=1992 => [
        #         "other_state_land", "mountain_reserves", "tea_development_forest",
        #         "forest", "private_forest_or_wasteland",
        #     ],
        #     water=[
        #         1992 => "lakes",
        #     ]
        # ),
        wlf = "$path/Data/Generated/Landcover/mus_wlf_shape.tif" =>
            ["cleared", "other"] => (; 
                native=2002 => "other",
                cleared=2002 => "cleared",
                abandoned=2002 => "other",
                urban=2002 => "other",
                forestry=2002 => "other",
                water=2002 => "other",
            ),
        homiisland = "$path/Data/Generated/Landcover/mus_landcover_2.tif" =>
            [
                "Continuous_urban",
                "Discontinuous_urban",
                "Forest",
                "Shrub_vegetation",
                "Herbaceaous_vegetation",
                "Mangrove",
                "Barren_Land",
                "Water",
                "Sugarcane",
                "Pasture",
                "",
                "Other_cropland",
            ] =>
            (;
                native=2017 => ["Forest", "Shrub_vegetation"],
                cleared=2017 => ["Sugarcane", "Pasture", "Other_cropland", "Herbaceaous_vegetation"],
                abandoned=2017 => ["Barren_Land", "Pasture", "Shrub_vegetation", "Forest"],
                urban=2017 => ["Continuous_urban", "Discontinuous_urban", "Herbaceaous_vegetation"],
                forestry=2017 => "Forest",
                water=2017 => "Water",
            ),
        # mascarine_birds_1 = "/home/raf/PhD/Mascarenes/maps/Mauritius/Studies_of_Mascarine_birds.tif" =>
            # ["Forestry"] => (; forestry=1984 => "Forestry",),
        # ),
        forest = "$path/Data/Selected/Mauritius/forest.tif" =>
            [
                "low",
                "medium",
                "high",
            ] => (;
                native=[
                    2020 => ["low", "medium", "high"],
                ],
                cleared=2020 => (&, :mask, (!, ["low", "medium", "high"])),
                abandoned=2020 => (&, :mask, (!, ["low", "medium", "high"])),
                urban=2020 => (&, :mask, (!, ["low", "medium", "high"])),
            ),
        forestry_1975 = "$path/maps/Mauritius/Studies_of_Mascarine_birds/31.png-1.png" =>
            (;
                cleared=1975 => "cleared_1975",
                forestry=1975 => "cleared_1975",
            ),
        vegetation_1975 = "$path/maps/Mauritius/Studies_of_Mascarine_birds/9.png-1.png" =>
            (;
                native=1975 => ["surviving_native", "mixed_native_and_plantation"],
                cleared=1975 => "cleared_1975",
                abandoned=1975 => "exotic_scrub",
                forestry=1975 => ["forest_plantation", "mixed_native_and_plantation"]
            ),
        forestry_1980 = "$path/maps/Mauritius/Studies_of_Mascarine_birds/49_2.png-1.png" =>
            (;
                native=[
                    1947 => ["native_1947_or_1980", "native_1980"],
                    1980 => ["native_1980"],
                ],
            ),
        forestry_1984 = "$path/maps/Mauritius/Studies_of_Mascarine_birds/48.png-1.png" =>
            (;
                cleared=1984 => "cleared_1973-1984",
                forestry=1984 => "cleared_1973-1984",
            ),
        # vegetation = "$path/Data/Selected/Mauritius/Undigitised/page33_mauritius_vegetation_colored.png" => (;),
        # fraser_1835_composite_etsy = "$path/Data/Selected/Mauritius/Undigitised/1835_fraser_composite_etsy.png" => (;
        #     cleared=(cleared_1835="cleared",), uncleared=(forest_1835="uncleared",),
        # ),
        # surveyor_general_1872_from_gleadow = "$path/Data/Selected/Mauritius/Undigitised/1872_surveyor_general_from_gleadow.jpg" => (;
        #     # cleared=(cleared_1850="cleared_1850-70", cleared_1870=("cleared_1850-70", "cleared_1870-72"), cleared_1872=("cleared_1850-70", "cleared_1870-72")),
        #     forest=[
        #         1850 => ["forest_1872", "cleared_1850-70", "cleared_1870-72"],
        #         1870 => ["forest_1872", "cleared_1870-72"],
        #         1872 => "forest_1872",
        #     ],
        # ),
    ), reu=(;
        # cadet_invasives="$path/Data/Selected/Reunion/Undigitised/cadet_invasives.jpg" => (;)),
        # atlas_vegetation = "$path/Data/Selected/Reunion/Undigitised/atlas_vegetation.jpg" => (;)),
        # atlas_ownership = "$path/Data/Selected/Reunion/Undigitised/atlas_ownership.jpg" => (;)),
        # atlas_1960_population = "$path/Data/Selected/Reunion/Undigitised/atlas_1960_population.jpg" => (;)),
        # "atlas_1960_agriculture" => "$path/Data/Selected/Reunion/Undigitised/atlas_agriculture_1960.jpg" => (;)),
        # atlas_population_1967 = "$path/Data/Selected/Reunion/Undigitised/atlas_population_1967.jpg" => (;)),
        atlas_1960_agriculture = "$path/Data/Selected/Reunion/Undigitised/atlas_agriculture_1960_2.jpg" => (;
            native=1960 => ["forest", "shrubland", "rock", "savannah", "geranium_discontinuous"],
            cleared=1960 => ["cane", "geranium_continuous", "tea", "geranium_discontinuous"],
            forestry=1960 => ["casuarina", "acacia", "cryptomeria", "labourdonassia"],
            urban=1960 => "urban",
            abandoned=1960 => ["forest", "shrubland", "rock", "savannah"],
            water=nothing,
        ),
        atlas_1815_agriculture = "$path/Data/Selected/Reunion/Undigitised/atlas_1815_agriculture.jpg" => (;
            native = 1815 => "forest",
            cleared = 1815 => ["geranium", "vanilla", "cane"],
            abandonned = 1815 => ["wasteland", "forest"]
        ),
        atlas_1780_agriculture = "$path/Data/Selected/Reunion/Undigitised/atlas_1780_agriculture.jpg" => (;
            native = 1780 => "native",
            cleared = 1780 => (!, "native"),
        ),
        # atlas_early_settlement = " => (;
        #     concessions=(conceded_1715_1765="conceded_1715-1765", conceded_1665_1715="concede_1665-1715"),
        # ),
    ))
end

lc = Raster("/home/raf/PhD/Mascarenes/Data/Generated/Landcover/mus_landcover.tif")
Rasters.rplot(lc)

_joinpath(path::String, filename) = joinpath(path, filename)
_joinpath(path::String, filename) = joinpath(path, filename)

function get_times(categories)
    times = Set{Int}()
    for cat in categories
        if cat isa Pair{Int}
            push!(times, cat[1])
        elseif cat isa Vector
            foreach(cat) do (time, val)
                push!(times, time)
            end
        end
    end
    return sort!(collect(times))
end

function as_namedtuples(categories::NamedTuple, category_names, times)
    filled_times = map(times) do time
        map(category_names) do k
            haskey(categories, k) || return nothing
            category = categories[k]
            if category isa Pair
                if last(category) isa Symbol
                    return nothing
                elseif first(category) == time
                    return last(category)
                else
                    return nothing
                end
            elseif isnothing(category)
                return nothing
            else
                i = findfirst(c -> first(c) == time, category)
                if isnothing(i)
                    return nothing
                else
                    return last(category[i])
                end
            end
        end
    end
end

function make_raster_slices(masks, categories;
    path="/home/raf/PhD/Mascarenes",
    files = define_map_files(; path),
    category_names=nothing,
)
    println("Generating raster slices...")
    # Copy duplicated file wrap points
    fn = joinpath(path, "Data/Selected/Mauritius/Undigitised/atlas_19C_land_use")
    cp(fn * ".csv", fn * "_2.csv"; force=true)

    # Load all rasters and make masks layers for all land-cover classes.
    # We get the numbers form the saved ".json" file for the project.
    rasters = map(files, masks[keys(files)]) do island_files, mask
        files = map(island_files) do (image_path, categories)
            raster_path = splitext(image_path)[1] * ".tif"
            raw = mask .* fix_order(Raster(raster_path))
            if categories isa Pair # Manually pass in the names
                original_names = categories[1]
                categories = categories[2]
            else
                json_path = splitext(image_path)[1] * ".json"
                data = JSON3.read(read(json_path), MapRasterization.MapSelection)
                original_names = data.settings.category_name
            end
            times = get_times(categories)
            string_timeline = as_namedtuples(categories, category_names, times)
            # Gapfill NamedTuples with nothing
            grouped = map(string_timeline) do categories
                _category_raster(raw, original_names, categories, mask)
            end
            (; raw, grouped, times, string_timeline, original_names)
        end
        alltimes = sort!(union(map(f -> f.times, files)...))
        @show alltimes
        timeline_dict = Dict{Int,Any}()
        for file in files
            for time in alltimes
                i = findfirst(==(time), file.times)
                if !isnothing(i)
                    if haskey(timeline_dict, time)
                        timeline_dict[time] = map(.|, timeline_dict[time], file.grouped[i])
                    else
                        timeline_dict[time] = file.grouped[i]
                    end
                end
            end
        end
        timeline_pairs = sort!(collect(pairs(timeline_dict)); by=first)
        stacks = map(identity, RasterStack.(last.(timeline_pairs)))
        if eltype(stacks) == Any
            return nothing
        else
            years = first.(timeline_pairs)
            time = Ti(years; sampling=Intervals(End()), span=Irregular(1500, last(years)))
            timeline = RasterSeries(stacks, time)
            return (; files, timeline)
        end
    end
end

function _category_raster(raster::Raster, layer_names::Vector, layers::NamedTuple, mask)
    layers = map(layers) do layer
        _category_raster(raster, layer_names, layer, mask)
    end
end
function _category_raster(raster::Raster, layer_names::Vector, layer_components::Vector, mask)::Raster{Bool}
    layers = map(l -> _category_raster(raster, layer_names, l, mask), layer_components)
    out = rebuild(Bool.(broadcast(|, layers...)); missingval=false) .& mask
    @assert missingval(out) == false
    return out
end
function _category_raster(raster::Raster, layer_names::Vector, xs::Tuple{<:Function,Vararg}, mask)::Raster{Bool}
    f, args... = xs
    vals = map(args) do layer
        _category_raster(raster, layer_names, layer, mask)
    end
    return map(f, vals...) .& mask
end
function _category_raster(raster::Raster, layer_names::Vector, category::Symbol, mask)::Raster{Bool}
    if category === :mask
        return mask
    else
        error(":$category directive not understood")
    end
end
function _category_raster(raster::Raster, layer_names::Vector, category::Nothing, mask)::Raster{Bool}
    return map(_ -> false, raster)
end
function _category_raster(raster::Raster, layer_names::Vector, category::String, mask)::Raster{Bool}
    I = findall(==(category), map(String, layer_names))
    if length(I) == 0
        error("could not find $category in $(layer_names)")
    end
    # Get all values matching the first category as a mask
    out = rebuild(Bool.(raster .== first(I)); missingval=false)
    # Add pixels for any subsequent categories
    foreach(I[2:end]) do i
        out .|= raster .== first(i)
    end
    @assert eltype(out) == Bool
    @assert missingval(out) == false
    return out .& mask
end
function _category_raster(raster::Raster, layer_names::Vector, x, mask)
    error("slice must be a NamedTuple, String or Vector{String}, got $x")
end

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
    # if save
    #     output = MapRasterization.MapSelection(cs)
    #     json_path = splitext(img_path)[1] * ".json"
    #     if isfile(json_path)
    #         backup_path = splitext(img_path)[1] * "_backup.json"
    #         cp(json_path, backup_path; force=true)
    #     end
    #     write(json_path, JSON3.write(output))
    # end
    return cs
end

load_image(img_path::String) = RGB{Float64}.(load(img_path) |> rotr90)

namedvector_raster(stack::RasterStack) = namedvector_raster(NamedTuple(stack))
function namedvector_raster(stack::NamedTuple{K}) where K
    map(stack...) do args...
        NamedVector{K}(args)
    end
end

function stripe_raster(ser::AbstractRasterSeries, states)
    map(ser) do raster
        stripe_raster(raster, states)
    end
end
function stripe_raster(raster::AbstractRaster, states; stripedim=X())
    A = rebuild(similar(raster, Int); missingval=0)
    for I in DimIndices(A)
        v = raster[I...]
        c = count(v)
        x = if c == 0
            0
        elseif c == 1
            states[findfirst(Tuple(v))]
        else
            # Stripe the possible values
            stripe = mod1(val(dims(I, stripedim)), c)
            findall(Tuple(v))[stripe]
        end
        A[I...] = x
    end
    return A
end

"""
    compile_timeline(transitions, series::RasterSeries{<:Raster{<:NamedVector,2}}; kw...)
    compile_timeline(transitions, raster::Raster{<:NamedVector,3}; kw...)

Compile a series of rasters into a coherent timeline,
following the logic defined in `transitions`.

Returns a `RasterStack` containing both the compiled time series, 
and the error associated with forcing `transitions`.

# Keywords

-`fill`: back fill empty pixels by taking `|` of 
    the previous and next non-empty pixels.
-`continuity`: assume continuity between single valued pixels
    accross time. If any intermediate pixels include that value
    as well as others, assume there were no changes. This may
    vary between a realistic assumtion and a very unrealistic
    assumtion, depending on categories and the length of time.
-`transitions`: force all transition logic, forwards and backwards,
    recording all transition errors to a separate raster.
"""
compile_timeline(transitions, ser::RasterSeries; kw...) =
    compile_timeline(transitions, Rasters.combine(ser, Ti); kw...)
function compile_timeline(transitions::NamedVector{K}, raster::Raster{T}; kw...
) where {K,T}
    states = NamedVector{K}(ntuple(identity, length(K)))
    # Get time on the first dimension for performance
    timeline = copy(raster)#permutedims(raster, (Ti(), X(), Y()))
    error = fill!(similar(timeline), zero(T))
    pixel_timeline_copy1 = fill!(timeline[X(1), Y(1)], zero(T))
    pixel_timeline_copy2 = fill!(timeline[X(1), Y(1)], zero(T))
    for xy in DimIndices(dims(raster, (X(), Y())))
        pixel_timeline = view(timeline, xy...)
        pixel_error = view(error, xy...)
        pixel_timeline_copy1 .= pixel_timeline
        pixel_timeline_copy2 .= pixel_timeline
        compile_timeline!(pixel_timeline, transitions; 
            pixel_error, pixel_timeline_copy1, pixel_timeline_copy2, kw...
        ) 
    end
    missingval = (timeline=zero(eltype(timeline)), error=zero(eltype(error)))
    x = RasterStack((; timeline, error); missingval)
    return rebuild(x; missingval)
end

function compile_timeline!(pixel_timeline, transitions; 
    indirect = indirect_logic(transitions), 
    reversed = reverse_transitions(transitions), 
    reversed_indirect = reverse_transitions(indirect),
    pixel_error = copy(pixel_timeline),
    pixel_timeline_copy1 = copy(pixel_timeline),
    pixel_timeline_copy2 = copy(pixel_timeline),
    fill=true, continuity=true, transition=true,
) 
    if fill
        fill_empty_times!(pixel_timeline)
    end
    if continuity
        remove_intermediate_uncertainty!(pixel_timeline, transitions, reversed, indirect, reversed_indirect)
    end
    if transition
        rev = lastindex(pixel_timeline):-1:1
        apply_transitions!(pixel_timeline_copy1, pixel_error, reversed, reversed_indirect)
        apply_transitions!(view(pixel_timeline_copy2, rev), view(pixel_error, rev), transitions, indirect)
        if continuity
            remove_intermediate_uncertainty!(pixel_timeline_copy1, transitions, reversed, indirect, reversed_indirect)
            remove_intermediate_uncertainty!(view(pixel_timeline_copy2, rev), reversed, transitions, reversed_indirect, indirect)
        end
        broadcast!(pixel_timeline, pixel_timeline_copy1, pixel_timeline_copy2) do a, b
            ands = map(&, a, b)
            if any(ands)
                ands
            else
                map(|, a, b)
            end
        end
        if continuity
            remove_intermediate_uncertainty!(pixel_timeline, transitions, reversed, indirect, reversed_indirect)
        end
    end
    return pixel_timeline, pixel_error
end

function fill_empty_times!(timeline::AbstractVector)
    last_non_empty_i = typemax(Int)
    last_non_empty_categories = zero(eltype(timeline))
    for i in eachindex(timeline)
        present_categories = timeline[i]
        any(present_categories) || continue # No category data for this time slice
        if last_non_empty_i < (i - 1)
            # There has been a gap in data, fill it with the combination of 
            # the last non empty categories and the present category
            fill_categories = map(|, last_non_empty_categories, present_categories)
            for n in last_non_empty_i+1:i-1
                timeline[n] = fill_categories 
            end
        end
        last_non_empty_i = i
        last_non_empty_categories = present_categories
    end
end

function remove_intermediate_uncertainty!(timeline, transitions, reversed, indirect, reversed_indirect)
    past = first(timeline)
    last_single_i = typemax(Int)
    last_single_categories = zeros(eltype(timeline))
    for i in eachindex(timeline)
        current_categories = timeline[i]
        # If we have got to a certain category, or the end
        (count(current_categories) == 1 || (i == lastindex(timeline))) || continue
        if last_single_i < (i - 1)
            # Check that all intermediate timesteps contain
            # the present category as one of the possibilities
            fillrange = last_single_i+1:i-1
            # Replace the intermediate uncertain categories
            for n in fillrange
                categories_at_n = timeline[n]
                shared = map(&, last_single_categories, categories_at_n, current_categories)
                if any(shared)
                    timeline[n] = shared
                else
                    cats  = _merge_all_possible(last_single_categories, categories_at_n, reversed)
                    final = _merge_all_possible(current_categories, cats, transitions)
                    if any(final)
                        timeline[n] = final
                    else
                        cats2  = _merge_all_possible(last_single_categories, categories_at_n, reversed_indirect)
                        final2 = _merge_all_possible(current_categories, cats2, indirect)
                        timeline[n] = final2
                    end
                end
            end
        end
        last_single_i = i
        last_single_categories = current_categories
    end
end

function _merge_all_possible(source, dest, transitions, indirect)
    direct = _merge_all_possible(source, dest, transitions)
    if any(direct)
        return direct
    else
        return _merge_all_possible(source, dest, indirect)
    end
end
function _merge_all_possible(source, dest, transitions)
    reduce(zip(source, transitions); init=zero(source)) do acc, (s, possible_dest)
        if s
            map(|, map(&, dest, possible_dest), acc)
        else
            acc
        end
    end
end

function apply_transitions!(timeline, errors, transitions, indirect)
    a = first(timeline)
    for i in firstindex(timeline)+1:lastindex(timeline)
        b = timeline[i]
        matching_b = _merge_all_possible(a, b, transitions, indirect)
        a, error = if any(matching_b)
            matching_b, zero(a)
        else
            if any(a)
                if any(b)
                    # Nothing matches and `a` is certain, so
                    # propagate `a` and return `b` in error
                    if count(a) == 1
                        a, b
                    else
                        # Never force-propagate uncertain states
                        map(|, a, b), map(|, a, b)
                    end
                else
                    # b is empty, just propagate `a`
                    a, zero(b)
                end
            else
                # `a` is empty, just propagate `b`
                b, zero(a)
            end
        end
        timeline[i] = a
        errors[i] = errors[i] .| error
    end
end


function indirect_logic(transitions::NamedVector{K}) where K
    states = NamedVector{K}(ntuple(identity, length(K)))
    return indirect_logic(states, transitions)
end
function indirect_logic(states, transitions)
    map(states) do s1
        map(states) do s2
            can_change(states, transitions, s1, s2)
        end
    end
end

function can_change(states, transitions, to, from, checked=(), path=())
    if logic[to][from]
        return true
    else
        return map(states) do s
            if s == to || s in checked
                false
            elseif transitions[to][s]
                can_change(states, transitions, s, from, (checked..., to), (path..., from))
            else
                false
            end
        end |> any
    end
end

function reverse_transitions(transitions::NamedVector{K}) where K
    map(NamedVector{K}(K)) do k
        map(transitions) do l
            l[k]
        end
    end
end
