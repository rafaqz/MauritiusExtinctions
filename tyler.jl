using Tyler, TileProviders
using GLMakie
using MakieDraw
tyler = Tyler.Map(Extents.extent(dems.rod); provider=Google(), resolution=(1000, 1500))
tyler.axis.aspect = 1.6
poly_canvas = GeometryCanvas{Polygon}(; figure=tyler.figure, axis=tyler.axis)
GeoJSON.write("/home/raf/PhD/Mascarenes/Data/Generated/RestrictedAreas/rod_cattle_walk.json", poly_canvas)
