using GeoInterface
using Colors, GeometryBasics
using GLMakie
using Makie
using CSV
using MapRasterization
using DataFrames
using GeoJSON
using Tables
includet("raster_common.jl")
includet("roads.jl")
includet("water.jl")

svg_path = "/home/raf/PhD/Mascarenes/Data/Selected/Mauritius/forest.svg"
json_path = splitext(svg_path)[1] * ".json"
if isfile(json_path)
    native_veg_poly = GeoJSON.read(read(json_path))
end
native_veg_rast = round.(Union{Int,Missing}, rasterize(native_veg_poly; to=dems.mus, fill=:category))
native_veg_mask = boolmask(native_veg_rast)
plot(native_veg_mask)
grade_fractions = (0.7, 0.5, 0.2)
native_density = mask(rebuild(map(native_veg_rast) do x
    ismissing(x) ? 0.0 : grade_fractions[x]
end; name=:native_density); with=dems.mus)
plot(native_density)

# using XML
# """
#     flatten(obj) -> Vector
#     flatten(obj, match::Pair) -> Vector
#     flatten(f::Function, obj) -> Vector

# Recursively flatten the children of an XML object based on the `match`
# which is a `Pair` of `Type => tag`, or a function `f` of a child object
# that returns a `Bool`.

# `flatten` with a single argument will simply flatten children by one level.
# """
# flatten(obj) = flatten(_ -> true, obj)
# function flatten(obj, type::Type{T}) where {T <: XML.AbstractXMLNode}
#     flatten(x -> x isa T, obj)
# end
# function flatten(obj, (typ, tag)::Pair)
#     flatten(x -> x isa typ && XML.tag(x) == tag, obj)
# end
# function flatten(f::Function, obj)
#     objs = map(children(obj)) do child
#         if f(child)
#             return [child]
#         else
#             return flatten(f, child)
#         end
#     end
#     return collect(Iterators.flatten(filter(!isempty, objs)))
# end

# function to_geoms(::Type{T}, p) where T
#     matrix = get_transform(p)
#     geoms = T[]
#     parts = split(p.d, " ")
#     local line = Point2{Float64}[]
#     s = nothing
#     while true
#         i = isnothing(s) ? iterate(parts) : iterate(parts, s)
#         if isnothing(i) # end of iteration 
#             if length(line) > 1
#                 push!(geoms, T(line))
#             end
#             break
#         end
#         x, s = i
#         if x == "M"
#             line = Point2{Float64}[]
#             a, s = iterate(parts, s)
#             b, s = iterate(parts, s)
#             push!(line, parsepoint(matrix, a, b))
#         elseif x == "L"
#             a, s = iterate(parts, s)
#             b, s = iterate(parts, s)
#             push!(line, parsepoint(matrix, a, b))
#         elseif x == "Z"
#             push!(line, first(line))
#             if length(line) > 1
#                 push!(geoms, T(line))
#             end
#         end
#     end
#     geoms
# end

# function parse_matrix(s)
#     m = parse.(Float64, split(s[8:end-1], ","))
#     [m[1] m[3] m[5];
#      m[2] m[4] m[6];
#      0.0  0.0  0.1 ;
#     ]
# end

# get_transform(p) = hasproperty(p, :transform) ? parse_matrix(p.transform) : nothing

# parsepoint(m::Nothing, a, b) = Point2(_parsefloats(a, b))
# parsepoint(m, a, b) = Point2(m * [_parsefloats(a, b)..., 1])
# _parsefloats(a, b) = [parse(Float64, a), parse(Float64, b)]

# # use = flatten(xml.root, Element => "use")
# # glyphs = map(children(first(children(first(children(xml.root)))))) do g
# #     if length(children(g)) > 0 
# #         path = first(children(g))
# #         occursin("C", path.d) ? missing : (g.id => to_geoms(Polygon, path))
# #     else
# #         missing
# #     end
# # end |> skipmissing |> collect |> Dict
# # use = filter(children(xml.root)) do e
# #     e isa Element && tag(e) == "g"
# # end
# # glyph_uses = map(first.(children.(use)), ) do u
# #     (href=getproperty(u, Symbol("xlink:href"))[2:end], pos=Point(parse(Float64, u.x), parse(Float64, u.y)))
# # end
# # pairs = map(glyph_uses) do gu
# #     gu => findfirst(x -> x.id == gu.href, glyphs)
# # end |> ps -> filter(p -> !isnothing(last(p)),  ps) 
# # href_polygons = map(pairs) do (u, n)
# #     Polygon(glyphs[n].geom .+ u.pos)
# # end

# xml = XML.Document(svg_path);

# allpaths = flatten(xml.root, Element => "path")

# color_strings = map(allpaths) do x
#     hasproperty(x, :fill) || return "none"
#     x.fill == "none" ? (hasproperty(x, :stroke) ? x.stroke : "none") : x.fill 
# end |> union

# category_paths = map(color_strings) do color_string
#     if color_string == "none"
#         color = RGB(0.0, 0.0, 0.0)
#     else
#         color = RGB(map(s -> parse(Float64, s) / 100, split(color_string[5:end-2], "%, "))...)
#     end
#     paths = filter(allpaths) do p
#         hasproperty(p, :fill) || return color_string == "none"
#         p.fill == "none" || return p.fill == color_string 
#         if hasproperty(p, :stroke) 
#             p.stroke == color_string
#         else
#             color_string == "none"
#         end
#     end
#     linestrings = Vector{Tuple{Float64,Float64}}[]
#     T = occursin("Z", first(paths).d) ? Polygon : LineString
#     geoms = map(paths) do p
#         if occursin("C", p.d)
#             return missing
#         else
#             to_geoms(T, p)
#         end
#     end |> x -> map(identity, Iterators.flatten(filter(!ismissing, x))) # Not collect, for better eltype
#     (; color, geoms)
# end;
# geoms = category_paths[1].geoms

# # CairoMakie.activate!()
# # colors = getproperty.(category_paths, :color)
# # fig = Figure(resolution = (1800,1600))
# # fig[1,1] = ax = Makie.Axis(fig)
# # ax.xreversed = false
# # ax.yreversed = true
# # display(fig)
# # Makie.heatmap(lookup(template, X), lookup(template, Y), template) 
# # c = category_paths[9]
# # (x -> x).(c.geoms)
# # lines!(ax, c.geoms; color=c.color)
# # selected = [3,6,7,8,9,10]
# # for i in selected
# #     c = category_paths[i]
# #     if first(c.geoms) isa Polygon
# #         poly!(ax, c.geoms; color=c.color)
# #     else
# #         lines!(ax, c.geoms; color=c.color)
# #     end
# # end
# # display(fig)
# # save("/home/raf/plot.pdf", fig, px_per_unit = 4)

# road_cats = [3, 9, 10]
# # template = broadcast(watermasks.mus, dems.mus) do w, d
# #     !ismissing(w) && w == 1 ? Base.nonmissingtype(eltype(dems.mus))(1000) : d
# # end
# template = broadcast(way_masks.mus[:ways_1880], dems.mus) do w, d
#     !ismissing(w) && w == 1 ? Base.nonmissingtype(eltype(dems.mus))(1000) : d
# end
# csv_path = splitext(svg_path)[1] * ".csv"
# pts = CSV.read(csv_path, DataFrame)
# geoms = reduce(vcat, getproperty.(category_paths[road_cats], :geoms))
# warp_points = MapRasterization.select_warp_points(geoms; template, guide=road_lines)
# warp_points = MapRasterization.select_warp_points(geoms; template, points=pts, guide=road_lines)
# warp_points
# # CSV.write(csv_path, warp_points)

# vegetation_cats = [6,7,8]
# warped_paths = map(enumerate(category_paths[vegetation_cats])) do (i, c)
#     geoms = MapRasterization.applywarp(c.geoms; template, points=pts)
#     (; color=c.color, category=i, geoms)
# end
# # Plot
# using GLMakie
# GLMakie.activate!()
# fig = Figure(resolution = (1800,1600))
# fig[1,1] = ax = Makie.Axis(fig)
# Makie.heatmap!(ax, lookup(template, X), lookup(template, Y), template; colormap=:magma) 
# for (w, c) in zip(warped_paths, category_paths[vegetation_cats])
#     Makie.poly!(ax, w.geoms; color=c.color)
# end

# warped_vegetation = map(enumerate(warped_paths)) do (i, c)
#     map(c.geoms) do g
#         (; geometry=g, category=i, color=c.color)
#     end
# end |> xs -> vcat(xs...)
# warped_vegetation
# GeoJSON.write(json_path, warped_vegetation)
# warped_vegetation = GeoJSON.read(read(json_path))

# # Raster
# raster_path = splitext(svg_path)[1] * ".tif"
# isfile(raster_path)
# isfile("/home/raf/PhD/Mascarenes/Data/Selected/Mauritius/forest.svg")
# vegetation_raster = zeros(Int, dims(template))
# foreach(warped_vegetation) do v
#     rasterize!(vegetation_raster, v; to=template, fill=:category, missingval=0, boundary=:touches)
# end
# Makie.heatmap(lookup(vegetation_raster, X), lookup(vegetation_raster, Y), parent(vegetation_raster[Band(1)]); colormap=:magma) 
# write(raster_path, vegetation_raster)
# vegetation_raster = Raster(raster_path)
# raster_path
