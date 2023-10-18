using Rasters, Shapefile, GLMakie, DataFrames, Tyler, Proj, GeoInterfaceMakie, GeoInterface, TileProviders, Extents, GeometryBasics
using GeoDataFrames, XLSX, CSV, GADM, Statistics, GeometryOps, TerminalPager
const GI = GeoInterface
using DelimitedFiles

hyde_lu = replace_missing(RasterStack("/home/raf/Data/HYDE/baseline/zip/2017AD_lu"))
hyde_pop = replace_missing(RasterStack("/home/raf/Data/HYDE/baseline/zip/2017AD_pop"))
Rasters.rplot(hyde_pop.uopp_2017AD)
human_lu = broadcast(+, hyde_lu.cropland2017AD, hyde_lu.rangeland2017AD, hyde_pop.uopp_2017AD)

island_summaries = CSV.read("/home/raf/Data/Extinction/GIF_island_summaries_9feb21.txt", DataFrame)
island_locations = Point.(zip(island_summaries.Long, island_summaries.Lat))
extracted_cyclone = extract(cyclone_rast, island_locations) |> DataFrame
extracted_lu = extract(human_lu, island_locations) |> DataFrame
extracted_lu |> pager

island_summaries.cyclone_count = extracted.count
island_summaries
sort!(island_summaries, :cyclone_count)
model = lm(@formula(cyclone_count ~ Extinct), island_summaries)
r2(model)
island_summaries |> pager

provider = Provider("https://storage.googleapis.com/earthenginepartners-hansen/tiles/gfc_v1.7/last_543/{z}/{x}/{y}.jpg")
Tyler.Map(Extent(X=(1, 10), Y=(1, 10)); provider)


world_borders_path = "/home/raf/Data/world-administrative-boundaries/world-administrative-boundaries.shp"
world_borders_tbl = Shapefile.Table(world_borders_path) |> DataFrame
world_borders = GeoInterface.convert.(Ref(GeometryBasics), world_borders_tbl.geometry);
poly!(world_borders)

# Extinct plants of the world
extinct_plant_path = "/home/raf/Data/Extinction/Plants/41559_2019_906_MOESM3_ESM.xlsx"
extinct_plants_xl = XLSX.readxlsx(extinct_plant_path)
extinct_plant_df = DataFrame(XLSX.eachtablerow(extinct_plants_xl["Supplementary_Data_1"]))
filter!(row -> row.List == "Extinct", extinct_plant_df)
localities = extinct_plant_df.Locality
locality_codes = Iterators.flatmap(s -> split(s, ","), skipmissing(localities)) |> Base.union |> sort
# locality_conversion = Dict{String,Any}()
# for lc in locality_codes
#     if !haskey(locality_conversion, lc)
#         println("Enter code for $lc")
#         clipboard(lc)
#         input = readline()
#         if !(input == "")
#             locality_conversion[lc] = input
#         end
#     end
# end

# locality_conversion["WAS"] = "USA, Washington"
# locality_table = DataFrame(
#     :source_name => first.(collect(pairs(locality_conversion))),
#     :iso_name => last.(collect(pairs(locality_conversion))),
# )

# Not Easter Island and Juan Fernand are in the same file

# locality_table.iso_name = rstrip.(locality_table.iso_name)
locality_table = CSV.read("extinct_plant_locality_table.csv", DataFrame)
sort!(locality_table, [:island, :iso_name])
CSV.write("extinct_plant_locality_table.csv", locality_table)

gadm_regions = map(locality_table.iso_name) do region
    map(split(region, "+")) do subregion
        args = String.(split(subregion, ", "))
        try
            border = GeoInterface.convert(GeometryBasics, GADM.get(args...).geom[1])
            println("Found $subregion")
            return border
        catch
            println("Could not find \"$subregion\"")
            return missing
        end
    end
end;
region_multipolygons = map(gadm_regions) do geoms
    if ismissing(geoms) || any(ismissing, geoms)
        missing
    else
        map(geoms) do mp
            if GI.trait(mp) isa GI.PolygonTrait 
                [mp]
            else
                collect(GI.getgeom(mp))
            end
        end |> Iterators.flatten |> collect |> GeometryBasics.MultiPolygon
    end
end;
extinction_counts = map(locality_table.source_name) do name
    count(x -> occursin(name, x), skipmissing(extinct_plant_df.Locality))
end
full_locality_table = copy(locality_table)
full_locality_table.extinctions = extinction_counts
full_locality_table.geometry = region_multipolygons;
cyclone_freq = zonal(mean, cyclone_rast; of=full_locality_table, boundary=:touches)
full_locality_table.cyclone_freq = cyclone_freq

simplified = map(full_locality_table.geometry) do geom
    GeometryOps.simplify(geom; tol=0.001)
end;
full_locality_table.simplified = simplified;
filter!(r -> !ismissing(r.geometry), full_locality_table)
Rasters.rplot(rebuild(cyclone_rast; missingval=0))
poly!(collect(skipmissing(full_locality_table.simplified)); 
    colormap=(:reds, 0.5), color=full_locality_table.extinctions
)
sort!(full_locality_table, :extinctions)
full_locality_table |> pager

island_table = filter(r -> r.island && r.extinctions > 0, full_locality_table)
sort!(island_table, :extinctions)
island_table |> pager

using GLM
model = lm(@formula(cyclone_freq ~ extinctions), island_table)

@profview GeometryOps.simplify(full_locality_table.geometry[1:10]; tol=0.001)
using GLMakie
figure = Figure()
axis = Axis(figure[1, 1]; xlabel="n extinctions", ylabel="n cyclones")
p = scatter!(axis, log.(island_table.extinctions), log.(island_table.cyclone_freq))

# locality_table.geometry = borders;
# missed = filter(row -> any(ismissing, row.geometry), locality_table)
# sort(missed, :iso_name)
# using GeoJSON
# Make GeoJSON write tables...

GeoInterfaceMakie.@enable GeoInterface.Wrappers.WrapperGeometry

# extinct_shape_path = "/home/raf/Data/Extinction/redlist_species_data_a76bbdff-c908-4507-9196-772c77e11c96/data_0.shp"
# extinct_shapes = Shapefile.Table(extinct_shape_path) |> DataFrame

# wsg = Proj.reproject.(extinct_shapes.geometry; source_crs=EPSG(4326), target_crs=EPSG(3857))
# raster = rasterize(count, wsg; res=10000, boundary=:touches);
# Rasters.rplot(parent(raster))
# tyler = Tyler.Map(Extents.Extent(; X=(1, 80), Y=(1, 40)); provider=Google(), scale=2);
# heatmap!(tyler.axis, parent(raster); transparency=true, colormap=(:inferno, 0.5))

# using ProfileView
# for i in 1:length(wsg)
    # GLMakie.poly!(tyler.axis, wsg[i]; color=(:red, 0.1), transparency=true)
# end

# extinct_rast = rasterize(count, extinct_shapes; res=1/24, boundary=:touches)
# extinct_rast = coverage(extinct_shapes; mode=sum, res=1/24, boundary=:touches)
# plot(extinct_rast .> 5)
# plot(extinct_rast)
# plot!(borders_shp; fillcolor=nothing)

using CSV, DataFrames, Statistics, Chain
cyclones_path = "/home/raf/Data/Cyclones/bq-results-20230411-153235-1681227446653.csv"
cyclones_df = @chain begin
    CSV.read(cyclones_path, DataFrame)
    sort!(_, [:sid, :iso_time])
end
wind = map(eachrow(cyclones_df)) do row
    wind = (
        row.bom_wind,
        row.cma_wind,
        row.ds824_wind,
        row.hko_wind,
        row.mlc_wind,
        row.nadi_wind,
        row.neumann_wind,
        row.newdelhi_wind,
        row.reunion_wind,
        row.td9635_wind,
        row.td9636_wind,
        row.tokyo_wind,
        row.usa_wind,
        row.wellington_wind,
        row.wmo_wind,
    )
    w = skipmissing(wind)
    any(!ismissing, w) ? mean(w) : NaN
end
cyclones_df.wind = wind
saffir_simpson_ts = 63
saffir_simpson_1 = 119
cyclones_df = filter(r -> r.wind >= saffir_simpson_ts, cyclones_df)

points = map(Point, cyclones_df.longitude, cyclones_df.latitude)
scatter!(p.axis, points; color=cyclones_df.wind, colormap=:thermal)
# webmercator_points = map(parent(GeometryOps.reproject(points; source_crs=EPSG(4326), target_crs=EPSG(3857)))) do p
#     Point2(GeoInterface.x(p), GeoInterface.y(p))
# end
# cyclones_df.webmercator_points = webmercator_points
cyclones_df.points = points

# intensity = wind ./ maximum(filter(!isnan, wind))
cyclone_groups = groupby(cyclones_df, :sid)
cyclone_paths = map(pairs(cyclone_groups)) do (key, subdf)
    if any(p -> p[1] < 0, subdf.points)
        if any(p -> p[1] >= 0, subdf.points)
            # Split the one line that crosses zero
            [LineString(map(p -> Point(p[1] + 360, p[2]), filter(p -> p[1] < 0, subdf.points))),
             LineString(filter(p -> p[1] >= 0, subdf.points))
            ]
        else # All under zero
            [LineString(map(p -> Point(p[1] + 360, p[2]), subdf.points))]
        end
    else
        [LineString(subdf.points)]
    end
end |> Iterators.flatten |> collect |> x -> map(identity, x)
# cyclone_paths_wm = map(pairs(cyclone_groups)) do (key, subdf)
#     if any(p -> p[1] < 0, subdf.webmercator_points)
#         if any(p -> p[1] >= 0, subdf.webmercator_points)
#             [
#                 LineString(filter(p -> p[1] >= 0, subdf.webmercator_points)),
#                 LineString(filter(p -> p[1] < 0, subdf.webmercator_points)),
#             ]
#         else# All under zero
#             [LineString(map(Point(p[1] + 360, p[2]), subdf.webmercator_points))]
#         end
#     else
#         [LineString(filter(p -> p[1] >= 0, subdf.webmercator_points))]
#     end
# end |> Iterators.flatten |> collect
Makie.lines(cyclone_paths)

cyclone_rast_s = rasterize(count, cyclone_paths; res=1, to=Extent(X=(0, 360), Y=(-90, 90)))
cyclone_rast_s = rebuild(cyclone_rast_s; missingval=nothing)
cyclone_rast = set(cat(cyclone_rast_s[X(181:360)], cyclone_rast_s[X(1:180)]; dims=X), X=>-180:179)
Rasters.rplot(rebuild(cyclone_rast; missingval=0))
# Makie.poly!(borders; stokewidth=1)

# Rasters.metadata(cyclone_rast)[:missed_geometries]
# Rasters.rplot(cyclone_rast)

frug_path = "/home/raf/Data/Extinction/Islands/Data_Occurrences_IslandFrugivores.txt"
frug_df = CSV.read(frug_path, DataFrame)
union(frug_df[!, "Island.ID"])
frug_island_path = "/home/raf/Data/Extinction/Islands/Data_Island_Characteristics.txt"
frug_island_df = CSV.read(frug_island_path, DataFrame)
frug_combined_df = innerjoin(frug_df, frug_island_df, on=Symbol("Island.ID"), matchmissing=:notequal, makeunique=true)
frug_combined_df.WeigeltID = map(x -> x == "NA" ? missing : parse.(Int, x), frug_combined_df[!, "Weigelt.ID"])

island_shapes_path = "/home/raf/Data/Extinction/Islands/GADM_islands_Weigelt_etal_17883/GADM_islands_Weigelt_etal_17883/GADM_islands_17883data.shp"
island_shapes_df = GeoDataFrames.read(island_shapes_path)
names(island_shapes_df)
island_data_path = "/home/raf/Data/Extinction/Islands/Weigelt_etal_2013_PNAS_islanddata.csv"
island_df = CSV.read(island_data_path, DataFrame)
island_combined_df = innerjoin(island_df, island_shapes_df, on=[:Archip => :ARCHIP])#, matchmissing=:equal)
names(island_df)
island_combined_df.cyclone_freq = zonal(mean, cyclone_rast; of=island_combined_df.geometry, boundary=:touches)
xs = DataFrames.subset(island_combined_df, 
    :Dist => ByRow(>=(300)),
    :Area => ByRow(x -> 25000 > x > 10),
)
shps = GeoInterface.convert.(Ref(GeometryBasics), combined_df.geometry);
Makie.poly(borders; stokewidth=1, color=:grey)
Makie.poly!(shps; color=:red)

combined_df = innerjoin(frug_combined_df, island_combined_df, on=[:WeigeltID => :ID], matchmissing=:notequal, makeunique=true)
extinct_df = combined_df[combined_df.Status .== "Extinct", :]
extinct_df.Island
extant_df = combined_df[combined_df.Status .== "Extant", :]
id_groups = groupby(combined_df, :WeigeltID)
id_group_kv = collect(pairs(id_groups))
group_ids = map(x -> x[1][1], id_group_kv)
group_subdfs = Dict(group_ids .=> last.(id_group_kv))
n_extinct = map(island_combined_df.ID) do id
    if haskey(group_subdfs, id)
        count(==("Extinct"), group_subdfs[id].Status)
    else
        missing
    end
end
island_combined_df.n_extinct = n_extinct
used_df = filter(r -> !ismissing(r.n_extinct), island_combined_df)
plotpoints = Point2.(filter(p -> !any(ismissing, p), collect(zip(used_df.cyclone_freq, used_df.n_extinct))))
fig = Figure()
ax = Axis(fig[1, 1]; xlabel="cyclone frequency", ylabel="extinctions")
Makie.scatter!(ax, plotpoints)
text = map(s -> filter(isascii, s), used_df.Island)
length(text)
Makie.text!(ax, plotpoints; text)

using GLM
model = lm(@formula(n_extinct ~ cyclone_freq), island_combined_df)
using StatsPlots
StatsPlots.plot(predict(model, (; cyclone_freq = 0:100)))

island_combined_df.cyclone_gre
    cyclones = zonal(mean, cyclone_rast; of=island_df.geometry, boundary=:touches) 

background = zonal(mean, cyclone_rast; of=island_combined_df.geometry, boundary=:touches) |> mean

p = lines!(ax, freqs_at_n.n_extinct_freq, freqs_at_n.n)
lines(freqs_at_n.n, freqs_at_n.under_n_extinct_freq)
# shps = GeoInterface.convert.(Ref(GeometryBasics), extinction_shapes);
# shps = GeoInterface.convert.(Ref(GeometryBasics), no_extinction_shapes);
# Makie.poly(shps)

mean(filter(!isnan, extinct_df.cyclone_freq))
mean(filter(!isnan, extant_df.cyclone_freq))
length(filter(!isnan, extinct_df.cyclone_freq))
length(filter(!isnan, extant_df.cyclone_freq))
plot(extant_df.cyclone_freq)
plot!(extinct_df.cyclone_freq)
maximum(filter(!isnan, extant_df.cyclone_freq))

Makie.poly(shps; 
    colorrange=(0, maximum(filter(!isnan, meanfreq))),
    colormap=:inferno,
    color=meanfreq,
)

tyler = Tyler.Map(Extents.Extent(; X=(-180, 180), Y=(-80, 80)); provider=Google(), scale=2);
using Colors
color = map(intensity) do i
    RGBA(1.0, 0.0, 0.0, i)
end
Makie.scatter!(tyler.axis, mp_web; 
    transparency=true, 
    color,
    colormap=:thermal,
    # colorrange=extrema(cyclones_df.season),
)
Makie.scatter!(tyler.axis, mp_web; 
    transparency=true, 
    color=wind,
    colormap=:thermal,
    colorrange=(20.0, maximum(filter(!isnan, wind))),
)
Makie.lines!(tyler.axis, cyclone_paths_wm)
