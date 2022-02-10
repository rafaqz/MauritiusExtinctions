@time using GLM
@time using Rasters
# using ProfileView, Profile
# using SnoopCompile
@time using Images
@time using Makie
@time using GLMakie
@time using Colors
using Rasters.LookupArrays
using Rasters: Band
# tinf = @snoopi_deep using Images
# fg = flamegraph(tinf)
# ProfileView.view(fg)

# img = load("/home/raf/PhD/Mauritius/Data/reunion_clearing.png")
# raster = Raster(rotr90(img), (Y, X))
# raster = aggregate(:center, raster, 4)
#
# elevpath = "/home/raf/PhD/Mauritius/Data/Norder/LS factor/DEM/DEM100x100_Resample.img"
# elevation = Raster(elevpath; missingval=-3.4028235f38)[Band(1), X(1:520)]
# elev = replace_missing(elevation, 0)
# elevrgb = RGB24.(replace_missing(elev ./ maximum(elev), 1))

# x = Makie.image(parent(raster));
# scene = x.plot.parent
# unknownpoints = []
# on(scene.events.mousebutton) do event 
#     x, y = Makie.mouseposition_px(scene)
#     if event.button == Mouse.left && event.action == Mouse.press
#         println("Click ", (x, y))
#         xo, yo = Observable(x), Observable(y)
#         # scatter!(xo, yo)
#         push!(unknownpoints, (x, y))
#     end
# end
# scene = Scene(camera=campixel!, raw=true)
# image!(scene, elevrgb)

Aunknown = reverse(replace_missing(soilraster, 0) ./ maximum(soilraster); dims=2)
x = Makie.heatmap(axes(Aunknown)..., parent(Aunknown));
scene = x.plot.parent
unknownpoints = []
on(scene.events.mousebutton) do event 
    x, y = Makie.mouseposition(scene)
    if event.button == Mouse.left && event.action == Mouse.press
        println("Click ", (x, y))
        xo, yo = Observable(x), Observable(y)
        # scatter!(xo, yo)
        push!(unknownpoints, (x, y))
    end
end
display(scene)

unknownpoints
unknowntable = (((x, y),) -> (x_unknown=Float64(x), y_unknown=Float64(y))).(unknownpoints)
CSV.write("unknown.csv", unknowntable)

Aknown = reverse(replace_missing(elevation) ./ maximum(elevation); dims=1)
x = Makie.heatmap(map(parent, reverse(dims(Aknown)))..., parent(Aknown));
scene = x.plot.parent
knownpoints = []
on(scene.events.mousebutton) do event 
    x, y = Makie.mouseposition(scene)
    if event.button == Mouse.left && event.action == Mouse.press
        println("Click ", (x, y))
        xo, yo = Observable(x), Observable(y)
        # scatter!(xo, yo)
        push!(knownpoints, (x, y))
    end
end



knownpoints
knowntable = (((x, y),) -> (x_known=Float64(x), y_known=Float64(y))).(knownpoints)
combined = merge.(knowntable, unknowntable)
CSV.write("combined_points.csv", combined)
x_model = lm(@formula(x_unknown ~ x_known), combined)
y_model = lm(@formula(y_unknown ~ y_known), combined)
xs = round.(Int, predict(x_model))
ys = round.(Int, predict(y_model))

xs = round.(Int, predict(x_model, knowntable))
ys = round.(Int, predict(y_model, knowntable))
setindex!(Aunknown, 1, xs[1], ys[1])
setindex!(Aunknown, 1, xs[2], ys[2])
setindex!(Aunknown, 1, xs[3], ys[3])

rpoints = ((x, y) -> (x_known = x, y_known=y)).(DimTable(Aknown).X, DimTable(Aknown).Y)
xs = round.(Int, predict(x_model, rpoints))
ys = round.(Int, predict(y_model, rpoints))

Afixed = Aknown .* 0
Afixed .= missing

for (Ik, Iu) in  zip(CartesianIndices(Aknown), zip(xs, ys))
    if checkbounds(Bool, Aunknown, Iu...)
        @show Iu
        Afixed[Ik] = Aunknown[Iu...]
    end
end

x = Makie.heatmap(axes(Afixed)..., parent(Afixed))

function fitrasters(raster; to)
    clicks1 = Tuple{Int,Int}[]
    clicks2 = Tuple{Int,Int}[]
    button
        scene = Scene(camera=campixel!, raw=true)
        display(scene)
        image!(scene, to)
        on(scene.events.mousebutton) do button
            x, y = scene.events.mouseposition_px[]
            println("Click ", x, y, scene.px_area[].widths)
            scene = Scene(camera=campixel!, raw=true)
            image!(scene, r1)
            xo, yo = Observable(x), Observable(y)
            scatter!(xo, yo)
            push!(clicks1, pos)
        end
    end
    @sync begin
        scene = Scene(camera=campixel!, raw=true)
        display(scene)
        image!(scene, raster)
        on(scene.events.mousebutton) do button
            pos = scene.events.mouseposition[]
            println("Click ", pos, scene.px_area[].widths)
            scene = Scene(camera=campixel!, raw=true)
            image!(scene, r1)
            # scatter!(map(i -> [i], pos)...)
            push!(clicks2, pos)
        end
    end
    return clicks1, clicks2
end

clicks1, clicks2 = fitrasters(raster; to=elevrgb)
