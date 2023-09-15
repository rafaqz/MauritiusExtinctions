tbl = Shapefile.Table("/home/raf/PhD/Mascarenes/Data/Claudia/Black River Gorges and other shapes/npcs.shp") |> DataFrame
tbl = Shapefile.Table("/home/raf/PhD/Mascarenes/Data/Claudia/Demo/GIS WILD LIFE FOUNDATION/SHAPE FILE/land_use_WGS_region.shp") |> DataFrame
union(tbl.LAND_USE)
swamps = filter(:LAND_USE => x -> !ismissing(x) && x == "Marsh or Swamp", tbl)
cleared = filter( :LAND_USE => 
  x -> !ismissing(x) && x in ("Other Plantation", "Sugar Cane", "Other PlanTation", "Tea Plantation", "Sugar cane"), tbl)
# swamps_rast = boolmask(swamps; to=masks.mus)
cleared_rast = boolmask(cleared; to=masks.mus)
other_rast = .!(cleared_rast)
lc = cleared_rast .* 1 .+ other_rast .* 2 .* masks.mus
write("/home/raf/PhD/Mascarenes/Data/Generated/Landcover/mus_wlf_shape.tif", lc)
