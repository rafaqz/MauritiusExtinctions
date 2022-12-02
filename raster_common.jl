using GeoJSON
using GADM
using Shapefile
using RasterDataSources
using Rasters
using Plots
using Rasters: Between, trim, Band
using Plots: plot, plot!

workdir = "/home/raf/PhD/Mascarenes"
datadir = joinpath(workdir, "Data")
outputdir = joinpath(datadir, "Generated")
distancedir = joinpath(outputdir, "Distances")

includet("functions.jl")
includet("lost_land_images.jl")

years = 1638, 1773, 1835, 1872, 1935, "present"
lc_years = 1638, 1773, 1835, 1872, 1935, "present"
lc_year_keys = map(y -> "lc_$y", lc_years)

island_keys = (; mus=:mus, reu=:reu)

gdal_borders = (
    mus=GADM.get("MUS").geom[1],
    reu=GADM.get("REU").geom[1],
    # rod=GADM.get("MUS").geom[1],
)
borders = (mus=GeoInterface.convert(MultiPolygon, gdal_borders.mus),
           reu=MultiPolygon([GeoInterface.convert(Polygon, gdal_borders.reu)]))
bbox = (
    # mus=((57.1, 57.9), (-20.6, -19.8)), # with islands
    mus=((57.1, 57.9), (-20.6, -19.949)),
    reu=((55.0, 56.0), (-22.0, -20.0)),
    # rod = (63.0, 64.0), (-20.0, -19.0),
)
tiles = getraster(SRTM; bounds=bbox.mus)
dem1 = Raster(tiles[1]; name=:DEM)
dem2 = Raster(tiles[2]; name=:DEM)
border_selectors = map(bbox) do bb
    mus=(X(Between(bb[1])), Y(Between(bb[2])), Band(1))
end
# Mauritius is right over the split in the tiles
m1 = view(dem1, border_selectors.mus...)
m2 = view(dem2, border_selectors.mus...)
mus_dem = replace_missing(trim(cat(m1, m2; dims=Y); pad=10))
# Plots.plot(mauritius_dem)
reu_tile  = getraster(SRTM; bounds=bbox.reu)[1]
reu_dem = replace_missing(read(trim(view(Raster(reu_tile), border_selectors.reu...); pad=10)))
# rod_tile  = getraster(SRTM; bounds=rod_bounds)[1]
# rod_dem = trim(view(dem3, border_selectors...); pad=10)
dems = (mus=mus_dem, reu=reu_dem)

mauritius_proj_dem = Raster("/home/raf/PhD/Mascarenes/Data/Norder/LS factor/DEM/DEM100x100_Resample.img"; crs=EPSG(3337))
mauritius_proj_dem = rebuild(mauritius_proj_dem; missingval=minimum(mauritius_proj_dem))
plot(mauritius_proj_dem)
