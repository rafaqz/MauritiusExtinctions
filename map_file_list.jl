const HOMI_CLASSES = [
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
    "Other_cropland"
]
const HOMI = HOMI_CLASSES => (;
    native=2017 => ["Forest", "Shrub_vegetation"],
    cleared=2017 => ["Barren_Land", "Sugarcane", "Pasture", "Other_cropland"],
    abandoned=2017 => ["Forest", "Shrub_vegetation", "Barren_Land", "Herbaceaous_vegetation"],
    urban=2017 => ["Continuous_urban", "Discontinuous_urban"],
    forestry=2017 => "Forest",
    water=2017 => "Water",
)

function _get_categories(image_path)
    json_path = splitext(image_path)[1] * ".json"
    data = JSON3.read(read(json_path), MapRasterization.MapSelection)
    original_names = data.settings.category_name
end

function define_map_files(; 
    path = "/home/raf/PhD/Mascarenes"
)
    # Here we define:
    # - all of the files we use
    # - what the land-cover categories are called for the file
    # - a map between those categories in our main land cover categories
    #
    # These use either a single String or Tuples of strings staring with a function.
    # The function will later be broadcasted over masks of the separate layers to combine them
    # Mostly is `|` which is "or" so we make a mask of values that are true in one or the other file
    file_details = (mus=(;
        atlas_dutch_period = "$path/Data/Selected/Mauritius/Undigitised/atlas_dutch_period.jpg" => _get_categories => (;
            native=[
                1600 => :mask,
                1709 => ["ebony_harvest", "undisturbed"],
                1721 => ["ebony_harvest", "undisturbed"],
            ],
            # undisturbed=1709 => "undisturbed",
            # disturbed=1709 => "ebony_harvest",
            cleared=1709 => "cleared",
            abandoned=1721 => "cleared",
        ),
        atlas_18C_land_use = "$path/Data/Selected/Mauritius/Undigitised/atlas_18C_land_use.jpg" => _get_categories => (;
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
        fraser_1835_from_gleadow = "$path/Data/Selected/Mauritius/Undigitised/1835_fraser_from_gleadow.jpg" => _get_categories => (;
            abandoned=1835 => "forest",
            native=1835 => "forest",
            cleared=1835 => (!, "forest"),
            # urban=1835 => (!, "forest"),
            # water=1835 => (!, "forest"),
        ),
        atlas_19C_land_use = "$path/Data/Selected/Mauritius/Undigitised/atlas_19C_land_use.jpg" => _get_categories => (;
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
        landcover_1965 = "$path/Data/Selected/Mauritius/Undigitised/mus_landuse_1965_100_hi_c.pdf" => _get_categories => (;
            native=1965 => ["Forest_natural", "Swamps", "Rock", "Scrub", "Savannah"],
            cleared=1965 => ["Sugar", "Vegetables", "Tea"],
            abandoned=1965 => ["Rock", "Scrub", "Savannah"],
            urban=1965 => "Built_up",
            forestry=1965 => "Forest_plantation",
            water=[
                1965 => "Reservoirs",
            ],
        ),
        # atlas_1992_vegetation = "$path/Data/Selected/Mauritius/Undigitised/atlas_1992_vegetation.jpg" => (;),
        atlas_1992_agriculture = "$path/Data/Selected/Mauritius/Undigitised/atlas_1992_agriculture.jpg" => _get_categories => (;
            native=1992 => "forest",
            cleared=1992 => [
                "cane", "forage", "tea", "market_gardens", "cane_or_market_gardens",
                "cane_or_tea", "tea_or_market_gardens", "cane_or_fruit", "pasture",
            ],
            urban=1992 => "urban",
            forestry=[
                1992 => "forestry",
            ],
            water=1992 => "lakes",
            abandoned=1992 => ["pasture", "forest"],
        ),
        # We need a second round with this file as the categories overlap
        # This will overwrite anything incorrect in the first file for the specified year
        atlas_19C_land_use_2 = "$path/Data/Selected/Mauritius/Undigitised/atlas_19C_land_use_2.jpg" => _get_categories => (;
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
        forestry_1975 = "$path/maps/Mauritius/Studies_of_Mascarine_birds/31.png-1.png" => _get_categories => (;
                cleared=1975 => "cleared_1975",
                forestry=1975 => "cleared_1975",
            ),
        vegetation_1975 = "$path/maps/Mauritius/Studies_of_Mascarine_birds/9.png-1.png" => _get_categories => (;
                native=1975 => ["surviving_native", "mixed_native_and_plantation"],
                cleared=1975 => "cleared_1975",
                abandoned=1975 => "exotic_scrub",
                forestry=1975 => ["forest_plantation", "mixed_native_and_plantation"]
            ),
        forestry_1980 = "$path/maps/Mauritius/Studies_of_Mascarine_birds/49_2.png-1.png" => _get_categories => (;
                native=[
                    1947 => ["native_1947_or_1980", "native_1980"],
                    1980 => ["native_1980"],
                ],
            ),
        forestry_1984 = "$path/maps/Mauritius/Studies_of_Mascarine_birds/48.png-1.png" => _get_categories => (;
                cleared=1984 => "cleared_1973-1984",
                forestry=1984 => "cleared_1973-1984",
            ),
        wlf = "$path/Data/Generated/Landcover/mus_wlf_shape.tif" => ["cleared", "other"] => (;
                cleared=2002 => "cleared",
            ),
        homiisland = "$path/Data/Generated/Landcover/mus_landcover.tif" => HOMI,
        forest = "$path/Data/Selected/Mauritius/forest.tif" => ["low", "medium", "high"] => (;
            native=1999 => ["low", "medium", "high"],#, 2020 => ["low", "medium", "high"]], # Added to fix any errors in HOMI, as this likely hansn't changed a lot
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
        atlas_early_settlement = "$path/Data/Selected/Reunion/Undigitised/atlas_early_settlement_cropped.jpg" => _get_categories => (;
            native=[1715 => ["forest", "concede_1665-1715", "conceded_1715-1765"], 
                    1765 => ["forest", "concede_1665-1715", "conceded_1715-1765"]],
            # abandoned=1765 => "concede_1665-1715", # Guess (maybe wrongly) that land was not abanondoned at 50-0 years after concession
            cleared=[1715 => "concede_1665-1715", 
                     1765 => ["concede_1665-1715", "conceded_1715-1765"]],
        ),
        atlas_1780_agriculture = "$path/Data/Selected/Reunion/Undigitised/atlas_1780_agriculture.jpg" => _get_categories => (;
            native = [
                1600 => :mask,
                1780 => ["native", "rock"]
            ],
            cleared = 1780 => (!, ["native", "rock"]),
        ),
        atlas_1815_agriculture = "$path/Data/Selected/Reunion/Undigitised/atlas_1815_agriculture.jpg" => _get_categories => (;
            native = 1815 => ["wasteland", "forest"],
            cleared = 1815 => ["wasteland", "geranium", "vanilla", "cane"],
            abandoned = 1815 => ["wasteland", "forest"],
        ),
        public_map ="/home/raf/PhD/Mascarenes/Data/Selected/Reunion/Undigitised/34.jpg" => _get_categories => (;
            native=1958=>["vegetation1", "vegetation2", "vegetation3"],
            cleared=1958=>["cleared1", "cleared2"],
            abandoned=1958=>["vegetation1", "vegetation2", "vegetation3"],
            urban=1958=>["cleared1", "cleared2"],
            forestry=1958=>["vegetation1", "vegetation2", "vegetation3"],
            water=1958=>["lake"],
        ),
        atlas_1960_agriculture = "$path/Data/Selected/Reunion/Undigitised/atlas_agriculture_1960_2.jpg" => _get_categories => (;
            native=1960 => ["forest", "shrubland", "rock", "savannah", "geranium_discontinuous"],
            cleared=1960 => ["shrubland", "rock", "savannah", "cane", "geranium_continuous", "tea", "geranium_discontinuous"],
            abandoned=1960 => ["forest", "shrubland", "rock", "savannah"],
            urban=1960 => "urban",
            forestry=1960 => ["casuarina", "acacia", "cryptomeria", "labourdonassia"],
        ),
        # atlas_1960_land_owndership = "$path/Data/Selected/Reunion/Undigitised/atlas_ownership.jpg" => (;
        #     native=1960 => "National_forest_office",
        #     forestry=1960 => "National_forest_office",
        #     urban=1960 => "Habitat",
        # ),
        # homiisland = "$path/Data/Generated/Landcover/reu_landcover.tif" => HOMI,
        # natpark = "$path/Data/Generated/NationalParks/reu.tiff" => ["national park"] => (;
        #     native=2021 => "national park"
        # )
        native = "$path/Data/Generated/reu_all_natives.tif" => ["native_remnant"] => (;
            native=2005 => "native_remnant",
        ),
        quantifying_invasion_landcover = "$path/Data/Selected/Reunion/Undigitised/aec13048-fig-0001-m.jpg" => _get_categories => (;
            native=2019 => "Extant_native_vegetation",
            abandoned=2019 => "Disturbed_secondary_vegetation",
            forestry=2019 => "Disturbed_secondary_vegetation",
            cleared=2019 => "Agricultural_areas",
            urban=2019 => "Artificial_areas",
        ),
    ), rod=(;
        # This was not at all accurate the urban areas may be decorative or suggestive
        # urban_1885 = "/home/raf/PhD/Mascarenes/Data/Selected/Rodrigues/rodrigues-mauritius-mascarene-islands-mascarenhas-archipelago-1885-GC0ATP.jpg" => _get_categories => (;
        #     urban=1885 => "urban",
        # ),
        gade_1 = "$path/Data/Selected/Rodrigues/Gade_1985_landcover.png" => ["native_remnant"] => (;
            native = [1700 => :mask, 1985 => :force => "native_remnant"],
        ),
        gade_2 = "$path/Data/Selected/Rodrigues/Gade_1985_landcover_2.png" => _get_categories => (;
            abandoned=1985 => ["reforested", "fallow_or_settled"],
            cleared=1985 => ["cultivated", "grazing", "fallow_or_settled"],
            urban=1985 => "fallow_or_settled", # Unclear how large the "fallow" part is
        ),
        # homiisland = "$path/Data/Generated/Landcover/rod_landcover.tif" => HOMI_CLASSES => (;
        #         native=2017 => ["Forest"],
        #         cleared=2017 => ["Barren_Land", "Pasture", "Shrub_vegetation", "Sugarcane", "Herbaceaous_vegetation", "Other_cropland"],
        #         abandoned=2017 => ["Barren_Land", "Pasture", "Shrub_vegetation", "Forest", "Herbaceaous_vegetation"],
        #         urban=2017 => ["Continuous_urban", "Discontinuous_urban"],
        #     ),
        rural_development_planning = "/home/raf/PhD/Mascarenes/Data/Selected/Rodrigues/cb5989en.png-06_colored.png" => _get_categories => (;
            native=2021 => ["nature_reserve", "forest", "riparian_vegetation"],
            abandoned=2021 => ["nature_reserve", "forest", "riparian_vegetation"], 
            cleared=2021 => ["silvo_pasture", "pasture", "agriculture", "agricultural_residential"],
            urban=[#= 2016 => "urban_residential", =#2021 => "urban_residential"], # force urban into homiisland
        ),
    ))
end
