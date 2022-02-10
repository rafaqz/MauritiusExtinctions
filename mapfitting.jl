@time using GLM
@time using Rasters
@time using Images
@time using Makie
@time using GLMakie
@time using Colors
using Rasters.LookupArrays
using Rasters: Band

function selectposclick!(fig, ax, selectedpos)
    on(events(ax.scene).mousebutton, priority = 0) do event
        if event.button == Mouse.left
            if event.action == Mouse.press
                pos = mouseposition(ax.scene)
                pos_px = Makie.mouseposition_px(fig.scene)
                if in(pos_px, ax.scene.px_area[])                    
                    ipos = round.(Int, pos)
                    x = findfirst(p -> p == ipos, selectedpos[])
                    if isnothing(x)
                        # Add a point
                        selectedpos[] = push!(selectedpos[], ipos)
                    else
                        # Remove a point
                        deleteat!(selectedpos[], x)
                    end
                    notify(selectedpos)
                end
            end
        end        
        return Consume(false)
    end
end

function selectmultiple(A)
    fig = Figure()
    ax = Makie.Axis(fig[1,1])
    selectedpos = Observable(Point2{Int}[])
    Makie.heatmap!(ax, A)
    Makie.scatter!(ax, selectedpos, color=:red)
    selectposclick!(fig, ax, selectedpos)
    on(selectedpos) do pos
        @show pos
    end
    display(fig)
    return selectedpos
end

Aunknown = reverse(replace_missing(soilraster, 0) ./ maximum(soilraster); dims=2)
Aknown = reverse(replace_missing(elevation) ./ maximum(elevation); dims=2)

knownpoints = selectmultiple(parent(Aknown))
unknownpoints = selectmultiple(parent(Aunknown))

unknowntable = (((x, y),) -> (x_unknown=Float64(x), y_unknown=Float64(y))).(unknownpoints)
knowntable = (((x, y),) -> (x_known=Float64(x), y_known=Float64(y))).(knownpoints) combined = merge.(knowntable, unknowntable)
CSV.write("known.csv", unknowntable)
CSV.write("unknown.csv", unknowntable)
CSV.write("combined_points.csv", combined)

knowntable = CSV.File("unknown.csv")


# Linear models
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
