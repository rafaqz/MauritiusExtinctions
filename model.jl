using DynamicGrids,
      Plots,
      Rasters,
      Setfield, 
      Shapefile,
      StaticArrays,
      StaticBitArrays

includet("rules.jl")

# Data input

elevpath = "/home/raf/PhD/Mauritius/Data/Norder/LS factor/DEM/DEM100x100_Resample.img"
elevation = Raster(elevpath; missingval=-3.4028235f38)[Band(1)]
lakes_path = "/home/raf/PhD/Mauritius/Data/Norder/LS factor/Lakes/lakes_all.shp"
Shapefile.Handle(lakes_path) |> plot
soiltypes = Shapefile.Handle("/home/raf/PhD/Mauritius/Data/Norder/K factor/SoilK.shp")
water = "/home/raf/PhD/Mauritius/Data/Water_Areas/Mauritius_Water_Areas.shp"
water_shp = Shapefile.Handle(water)
plot(elevation)

elevation = Raster(elevpath; mappedcrs=EPSG(4326), missingval=-3.4028235f38)[Band(1)]
res = resample(elevation; crs=EPSG(4326))
plot(res)
p = plot(elevation; ylims=(-25, -10), xlims=(50, 70))
p = plot(elevation)
plot!(p, mauritius_border)
plot(Shapefile.Handle(rivers_path))
plot(elevation)
plot(soiltypes)

landuse_dir = "/home/raf/PhD/Mauritius/Data/Norder/C factor/"
years = 1638, 1773, 1835, 1872, 1935, "present"
landuse_shapefiles = map(years) do year
    path = joinpath(landuse_dir, string(year, ".shp"))
    Shapefile.Handle(path) 
en
landuse_snapshots = map(landuse_shapefiles, years) do shf, y
    landuse = zeros(Int, dims(elevation))
    # The forested/cleared order swaps after the first three files
    shapes = y in (1638, 1773, 1835) ? reverse(shf.shapes) : shf.shapes
    for (n, shape) in enumerate(shapes)
        rasterize!(landuse, shape; fill=n)
    end
    landuse
end

landuse_plots = map(landuse_snapshots, years) do sn, y
    plot(sn; title=string(y))
end
plot(landuse_plots...; size=(1400,1000))
p = plot(elevation)
plot(landuse_plots[2:3]...; opacity=0.3)

savefig("early_clearing.png")

land_use_category = (;
    forrested=1,
    cleared=2,
    urban=3,
)

# Grid initialisation
S = 100, 100
species = (SBitArray(UInt16, rand(Bool, 4, 4)) for i in 1:S[1], j in 1:S[2])
interaction_matrix = SArray(rand(10))
lu_response_matrix = map(_ -> (; forrested=rand(), cleared=rand(), urban=rand()), interaction_matrix)
hunting_susceptibility = SArray(rand(size(interaction_matrix))
landuse_susceptibility = SArray(rand(size(interaction_matrix))
landuse = zeros(Bool, S)

land_use_effect = Cell{Tuple{:S,:LU},:S}() do data, (p, l), I
    map(interaction_matrix, s) do i, s
        if i == zero(i)
            s
        else
            rand() < i ? false : s
        end
    end
end

species_interaction_effect = Cell{:S}() do data, s, I
    map(interaction_matrix, s) do i, s
        if i == zero(i)
            s
        else
            rand() < i ? false : s
        end
    end
end

rules = luc, land_clearing_effect, species_interaction_effect

output = ArrayOutput((species=species, land_use=land_use))

sim!(rules)

