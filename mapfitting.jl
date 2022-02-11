@time using GLM
@time using Rasters
@time using Images
@time using Makie
@time using GLMakie
@time using Colors
using Rasters.LookupArrays
using Rasters: Band

function manualwarp(As::Raster...; to::Raster)
    As = map(A -> reorder(A, ForwardOrdered), As)
    A1 = first(As)
    to = reorder(to, ForwardOrdered)
    # map(A -> size(A) == size(first(As)), As) || throw(ArgumentError("Intput raster sizes are not the same"))
    fig = Figure()
    ax1 = Makie.Axis(fig[1,1]; title="Source raster with known crs/resolution - `to` kw")
    ax2 = Makie.Axis(fig[1,2]; title="First raster with unknown crs/resolution")
    knownpoints = selectmultiple(parent(to), fig, ax1)
    unknownpoints = selectmultiple(parent(A1), fig, ax2)
    screen = display(fig)
    println("Select points in rasters, then close the window")
    while screen.window_open[] 
        sleep(0.1)
    end
    models = _fitlinearmodels(A1, to, (knownpoints, unknownpoints))
    warped = map(A -> linearwarp(A; to, models), As)
    display(Makie.heatmap(map(parent, dims(first(warped)))..., parent(first(warped))))
    if length(warped) == 1
        return first(warped)
    else
        return warped
    end
end

function selectposclick!(fig, ax, sct, positions)
    dragging = Ref(false)
    idx = Ref(0)
    on(events(fig).mousebutton, priority = 2) do event
        pos = Makie.mouseposition(ax.scene)
        ipos = round.(Int, pos)
        pos_px = Makie.mouseposition_px(fig.scene)
        if event.button == Mouse.left
            if event.action == Mouse.press
                if pos_px in ax.scene.px_area[]                    
                    plt, i = pick(fig.scene, pos...)
                    if plt == sct
                        @show i pos
                        idx[] = i
                    else
                        @show pos
                        push!(positions[], ipos)
                        idx[] = lastindex(positions[])
                        notify(positions)
                    end
                    dragging = true 
                end
            elseif event.action == Mouse.release
                dragging = false
            end
        elseif event.button == Mouse.right
            if pos_px in ax.scene.px_area[]                    
                x = findfirst(p -> p == ipos, positions[])
                if !isnothing(x)
                    # Remove a point
                    deleteat!(positions[], x)
                end
                notify(positions)
            end
        end
        return Consume(dragging[])
    end
    on(events(fig).mouseposition, priority = 2) do mp
        pos = Makie.mouseposition(ax.scene)
        ipos = round.(Int, pos)
        if dragging[]
            positions[][idx[]] = ipos
            notify(positions)
            return Consume(true)
        end
        return Consume(false)
    end
end

function selectmultiple(A, fig, ax)
    dragging = Ref(false)
    Makie.heatmap!(ax, A)
    positions = Observable(Point2{Float32}[])
    sct = Makie.scatter!(ax, positions, color=1:30, colormap=:reds)
    selectposclick!(fig, ax, sct, positions)
    return positions
end

function _fitlinearmodels(A, to, (knownpoints, unknownpoints))
    length(knownpoints[]) == length(unknownpoints[]) || error("Number of selected points must be the same for each raster")
    knowntable = (((x, y),) -> (x_known=Float64(x), y_known=Float64(y))).(knownpoints[])
    unknowntable = (((x, y),) -> (x_unknown=Float64(x), y_unknown=Float64(y))).(unknownpoints[])
    combined = merge.(knowntable, unknowntable)
    x_model = lm(@formula(x_unknown ~ x_known), combined)
    y_model = lm(@formula(y_unknown ~ y_known), combined)
    return x_model, y_model
end

function linearwarp(A; to, points=nothing, models::Union{Nothing,Tuple}=nothing)
    x_model, y_model = if isnothing(models) 
        isnothing(points) && error("pass either `points::Tuple` to fit or fitted `models::Tuple`")
        _fitlinearmodels(A, to, points)
    else
        models
    end

    knownpoints = vec(collect((x_known = x, y_known=y) for (x, y) in Tuple.(CartesianIndices(to))))
    xs = round.(Int, predict(x_model, knownpoints))
    ys = round.(Int, predict(y_model, knownpoints))
    Awarped = Raster(similar(to, Union{Missing,eltype(A)}); missingval=missing)
    Awarped .= missingval(Awarped)
    for (Ik, Iu) in  zip(CartesianIndices(to), CartesianIndex.(zip(xs, ys)))
        if checkbounds(Bool, A, Iu)
            Awarped[Ik] = A[Iu]
        end
    end
    return Awarped
end
