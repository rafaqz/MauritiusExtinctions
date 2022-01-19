using Rasters
using GeoInterface
using Shapefile
using Plots

elevpath = "/home/raf/PhD/Mauritius/Data/Norder/LS factor/DEM/DEM100x100_Resample.img"
lakespath = "/home/raf/PhD/Mauritius/Data/Norder/LS factor/Lakes/lakes_all.shp"
soiltypespath = "/home/raf/PhD/Mauritius/Data/Norder/K factor/SoilK.shp"
rainfallpath = "/home/raf/PhD/Mauritius/Data/Norder/R factor/r_annual.img"
landusedir = "/home/raf/PhD/Mauritius/Data/Norder/C factor/"

elevationraster = Raster(elevpath; missingval=-3.4028235f38)[Band(1)]
rainfallraster = Raster(rainfallpath)

lakes = Shapefile.Handle(lakes_path)
lakesraster = Raster{Int}(undef, dims(r); missingval=0)
for i in eachindex(lakes.shapes)[1:end-1]
    rasterize!(lakesraster, lakes.shapes[i]; fill=i)
end

soiltypes = Shapefile.Handle(soiltypespath)
soilraster = lakesraster .* 0
for i in eachindex(soiltypes.shapes)
    rasterize!(soilraster, soiltypes.shapes[i]; fill=i)
end

plot(soilraster)


# Rasterize waterways

years = 1638, 1773, 1835, 1872, 1935, "present"
landuse_shapefiles = map(years) do year
    path = joinpath(landuse_dir, string(year, ".shp"))
    Shapefile.Handle(path) 
end

landuse_snapshots = map(landuse_shapefiles, years) do shapefile, year
    landuse = zeros(Int, dims(elevation))
    # The forested/cleared order swaps after the first three files
    shapes = year in (1773, 1835) ? shapefile.shapes : reverse(shapefile.shapes)
    for (n, shape) in enumerate(shapes)
        rasterize!(landuse, shape; fill=n)
    end
    return landuse
end
plot(plot.(landuse_snapshots; clims=(0, 2), c=:viridis)...; size=(2000,2000))

#=
Projection    LAMBERT
Units         METERS
Zunits        NO
Xshift        0.0
Yshift        0.0
Parameters    6378137.0  6356752.314245179
 -20  0  0.0 /* 1st standard parallel
  57  0  0.0 /* 2nd standard parallel
  57 31 19.552232 /* central meridian
 -20 11 42.156258 /* latitude of projection's origin
1000000.0 /* false easting (meters)
1000000.0 /* false northing (meters)
=#
# 57.52209784
# -20.19504340

proj = WellKnownText{GeoFormatTypes.CRS, String}(
    GeoFormatTypes.CRS(), 
"""
PROJCS[\"World_Sinusoidal\",
    GEOGCS[\"WGS 84\",
        DATUM[\"WGS_1984\",
            SPHEROID[\"WGS 84\",6378137,298.257223563, 
                AUTHORITY[\"EPSG\",\"7030\"]
            ], 
            AUTHORITY[\"EPSG\",\"6326\"]
        ],
        PRIMEM[\"Greenwich\",0],
        UNIT[\"Degree\",0.0174532925199433]
    ], 
    PROJECTION[\"Sinusoidal\"],
    PARAMETER[\"central_meridian\",57.52209784],
    PARAMETER[\"latitude_of_origin\",20],
    PARAMETER[\"false_easting\",1000000.0],
    PARAMETER[\"false_northing\",3240000.0],
    UNIT[\"metre\",1,
        AUTHORITY[\"EPSG\",\"9001\"]
    ],
    AXIS[\"Easting\",EAST],
    AXIS[\"Northing\",NORTH]
]")
"""
)
