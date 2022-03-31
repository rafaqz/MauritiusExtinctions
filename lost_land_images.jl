lostland_image_classes = (
    mus=(
        phase=(
            magenta = "deforested before 1807",
            dark_red = "deforested 1807-1835",
            blue = "deforested 1835-1910",
            dark_blue = "deforested 1910-1947",
            green = "deforested 1947-1970",
            red = "deforested since 1970",
            yellow = "scrub with native remnants",
            dark_green = "native vegetation",
            cyan = "swamp",
        ),
        rem=(
            red = "native at least 20%",
            green = "native 50%+",
            blue = "reservoirs",
            cyan = "cleared",
        ),
        veg=(
            red = "semi-dry evergreen forest",
            green = "open dry palm-rich woodland",
            blue = "wet forest",
            cyan = "Pandanus swamp",
            magenta = "mossy rainforest",
            yellow = "mangrove",
            dark_red = "wetland vegetation",
        ),
    ),
    reu=(
        phase=(
            blue = "C17",
            red = "C18",
            green = "C19",
            yellow = "C20",
            cyan = "not cleared/colonised",
        ),
        rem=(
            red = "Lowland mixed evergreen semi-dry forest",
            green = "Lowland mixed evergreen tropical forest",
            blue = "Transitional mixed evergreen tropical rainforest",
            yellow = "Montane mixed evergreen subtropical rainforest",
            magenta = "Tree-heather formations",
            cyan = "Acacia-dominated montane rainforest",
            dark_red = "Very wet screw-pine formations",
            dark_green = "Upper montane temperate heaths",
            dark_blue = "Wetland vegetation",
            dark_yellow = "Largely unvegetated recent lava flows",
            dark_magenta = "Cleared",
        ),
        veg=(
            dark_cyan = "Palm-bejoin savanna",
            magenta = "Palm-rich dry forest",
            red = "Lowland mixed evergreen semi-dry tropical forest",
            yellow = "Lowland mixed evergreen tropical rainforest", 
            cyan = "Transitional mixed evergreen tropical rainforest",
            green = "Montane mixed evergreen subtropical rainforest",
            dark_blue = "Tree-heather formations",
            dark_green = "Acacia-dominated montane rainforest",
            dark_red = "Very wet screw-pine formations",
            blue = "Upper montane temperate heaths",
            dark_yellow = "Wetland vegetation",
            dark_magenta = "Unvegetated recent lava flows",
        ),
        # bulbul=load(joinpath(datadir, "LostLand/Maps/page265_reunion_bulbul_colored.png")),
    )
)

# Cleaning re-colored maps
lostland_image_paths = (
    mus=(;
        phase=joinpath(datadir, "LostLand/Maps/page132_mauritius_phases_colored.png"),
        rem=joinpath(datadir, "LostLand/Maps/page159_mauritius_remnants_colored.png"),
        veg=joinpath(datadir, "LostLand/Maps/page33_mauritius_vegetation_colored.png"),
        # kestrel=load(joinpath(datadir, "LostLand/Maps/page252_mauritius_kestrel_colored.png")),
    ),
    reu=(;
        phase=joinpath(datadir, "LostLand/Maps/page145_reunion_phases_colored.png"),
        rem=joinpath(datadir, "LostLand/Maps/page183_reunion_remnants_colored.png"),
        veg=joinpath(datadir, "LostLand/Maps/page36_reunion_vegetation_colored.png"),
        # bulbul=load(joinpath(datadir, "LostLand/Maps/page265_reunion_bulbul_colored.png")),
    )
)
