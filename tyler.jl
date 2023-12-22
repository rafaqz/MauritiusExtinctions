using Tyler, TileProviders
using GLMakie
using MakieDraw
tyler = Tyler.Map(Extents.extent(dems.rod); provider=Google(), resolution=(1000, 1500))
tyler.axis.aspect = 1.6
poly_canvas = GeometryCanvas{Polygon}(; figure=tyler.figure, axis=tyler.axis)
GeoJSON.write("/home/raf/PhD/Mascarenes/Data/Generated/RestrictedAreas/rod_cattle_walk.json", poly_canvas)


peiades_shp1 = "/home/raf/PhD/Mascarenes/Data/Peiades/Classif_2018_Peiades_Code1/Classif_2018_PL_final_Code1_communes.shp"
peiades_shp2 = "/home/raf/PhD/Mascarenes/Data/Peiades/Classif_2018_Peiades_Code2/Classif_2018_PL_final_Code2_communes.shp"
peiades_shp3 = "/home/raf/PhD/Mascarenes/Data/Peiades/Classif_2018_Peiades_Code3/Classif_2018_PL_final_Code3_communes.shp"
source_crs = GeoInterface.crs(Shapefile.Table(peiades_shp1))
p1 = Shapefile.Table(peiades_shp1) |> DataFrame
p2 = Shapefile.Table(peiades_shp2) |> DataFrame
p3 = Shapefile.Table(peiades_shp3) |> DataFrame
using GeometryOps
union(p1.Niveau1)
labels = map(union(p3.Niveau3)) do s
    in, out = Pipe(), Pipe()
    run(pipeline(`trans -brief $s`), in, out)
    close(in.in); close(out.in)
    s => String(read(out))
end
labels
union(p3.Niveau3)
forestry = filter(r -> ismissing(r.Niveau3) ? false : r.Niveau3 == "Plantation forestiere", p3)
forestry_webmerc = GeometryOps.reproject(forestry.geometry; source_crs, target_crs=EPSG(3857))
forestry_rast = boolmask(forestry_webmerc; res=30, verbose=false, boundary=:touches)

p = Makie.plot(forestry_rast; colormap=:magma)
Makie.plot!(p.axis, borders.reu)
save("pleiades_forestrry.png", p)
using Tyler, TileProviders, MapTiles
using GeoInterfaceMakie
ext = reduce(Extents.union, map(GeoInterface.extent, forestry_latlon))

tyler = Tyler.Map(ext, MapTiles.WebMercator(); provider=Google(:satelite))
Makie.image!(tyler.axis, forestry_rast; transparency=true, colormap=(:reds, 0.5))


rs = rebuild(resample(striped_raw.rod[Ti=At(2021)]; crs=EPSG(3857)); missingval=0)
slope.rod
Makie.heatmap!(tyler.axis, rebuild(resample(rs; crs=EPSG(3857), size=(1000, 1000)), missingval=0), colormap=(:magma, 0.5), transparency=true, opacity=0.2)
rs = rebuild(resample(slices.rod.files.gade_1.raw; crs=EPSG(3857), size=(1000, 1000)) .== 1; missingval=0)
sl = rebuild(resample(replace_missing((slices.rod.files.gade_1.raw .== 1) .* slope_stacks.rod.slope, 0); crs=EPSG(3857), size=(1000, 1000)); missingval=0)
sl = rebuild(resample(replace_missing(slope_stacks.rod.slope, 0); crs=EPSG(3857), size=(1000, 1000)); missingval=0)
Makie.heatmap!(tyler.axis, sl)
Makie.heatmap!(tyler.axis, rs; colormap=(:magma, 0.5), transparency=true, opacity=0.2)
