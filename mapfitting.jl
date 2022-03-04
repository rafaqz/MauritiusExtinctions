module RasterUtils

using Rasters
using GLM
using Makie
using GeometryBasics
using Colors

using Makie
using Rasters.LookupArrays

function manualwarp(As...; template::Raster, points=nothing, missingval=missing)
    A1 = first(As)
    if A1 isa Raster
        As = map(A -> reorder(A, ForwardOrdered), As)
    end
    template = reorder(template, ForwardOrdered)
    if isnothing(points)
        points = select_common_points(A1; template)
    end
    warped = _warp_from_points(As, template, points, missingval)
    # Show updated heatmap
    # display(Makie.heatmap(map(parent, dims(first(warped)))..., parent(first(warped))))
    # display(Makie.image(parent(first(warped))))
    if length(warped) == 1
        return first(warped)
    else
        return warped
    end
end

function _warp_from_points(As::Tuple, template, points, missingval)
    models = _fitlinearmodels(points)
    return map(A -> linearwarp(A; template, models, missingval), As)
end

select_common_points(A; template) = _select_common_points(A, template) 

_select_common_points(A, template::Raster) = _select_common_points(A, parent(reorder(template, ForwardOrdered)))
function _select_common_points(A, template)
    # map(A -> size(A) == size(first(As)), As) || throw(ArgumentError("Intput raster sizes are not the same"))
    fig = Figure()
    ax1 = Makie.Axis(fig[1,1]; title="Source raster with known crs/resolution - `template` kw")
    ax2 = Makie.Axis(fig[1,2]; title="First raster with unknown crs/resolution")
    knownpoints = selectmultiple(parent(template), fig, ax1)
    unknownpoints = selectmultiple(A, fig, ax2)
    screen = display(fig)
    println("Select points in rasters, then close the window")
    while screen.window_open[] 
        sleep(0.1)
    end
    length(knownpoints[]) == length(unknownpoints[]) || error("Number of selected points must be the same for each raster")
    knowntable = (((x, y),) -> (x_known=Float64(x), y_known=Float64(y))).(knownpoints[])
    unknowntable = (((x, y),) -> (x_unknown=Float64(x), y_unknown=Float64(y))).(unknownpoints[])
    return merge.(knowntable, unknowntable)
end

function manualinput(A::Raster; points=Point2{Float32}[])
    points = Point2{Float32}.(points)
    A = reorder(A, ForwardOrdered)
    fig = Figure()
    ax = Makie.Axis(fig[1,1]; title="Source image")
    ax.aspect = AxisAspect(1)
    Makie.heatmap!(ax, parent(A))
    positions = Observable(points)
    sct = Makie.lines!(ax, positions; color=:red)
    sct = Makie.scatter!(ax, positions; color=:red)
    dragselect!(fig, ax, sct, positions, size(A); caninsert=true)
    screen = display(fig)
    println("Select polygons in rasters, then close the window")
    while screen.window_open[] 
        sleep(0.1)
    end
    return positions[]
end

function dragselect!(fig, ax, sct, positions, pixelsize; caninsert=false)
    # Get pixel click accuracy from the size of visable heatmap.
    accuracy = lift(ax.finallimits) do fl
        round(Int, maximum(fl.widths) / 200)
    end
    dragging = Ref(false)
    idx = Ref(0)
    # Mouse down event
    on(events(fig).mousebutton, priority = 2) do event
        pos = Makie.mouseposition(ax.scene)
        ipos = round.(Int, pos)
        pos_px = Makie.mouseposition_px(fig.scene)
        # Add points with left click
        if event.button == Mouse.left
            if event.action == Mouse.press
                if pos_px in ax.scene.px_area[]                    
                    plt, i = pick(fig.scene, pos...)
                    idx[] = i
                    insert = false
                    found = pointnear(positions[], ipos, accuracy[]) do i
                        if isnothing(i) 
                            return nothing
                        else
                            idx[] = i
                            true
                        end
                    end
                    if isnothing(found)
                        if caninsert && length(positions[]) > 1
                            lastp = positions[][end]
                            # Search backwards so we preference recent lines
                            for i in eachindex(positions[])[end-1:-1:1]
                                p = positions[][i]
                                online = ison(Line(Point(lastp...), Point(p...)), Point(ipos...), accuracy[])
                                @show online
                                if online
                                    insert = true
                                    idx[] = i + 1
                                    insert!(positions[], i + 1, ipos)
                                    break
                                end
                                lastp = p
                            end
                        end
                        if !insert
                            push!(positions[], ipos)
                            idx[] = lastindex(positions[])
                        end
                    end
                    notify(positions)
                    dragging = true 
                end
            elseif event.action == Mouse.release
                dragging = false
            end
        # Delete points with right click
        elseif event.button == Mouse.right
            if pos_px in ax.scene.px_area[]                    
                pointnear(positions[], ipos, accuracy[]) do i
                    isnothing(i) || deleteat!(positions[], i)
                end
                notify(positions)
            end
        end
        return Consume(dragging[])
    end
    # Mouse drag event
    on(events(fig).mouseposition, priority = 2) do mp
        if dragging[]
            pos = Makie.mouseposition(ax.scene)
            ipos = round.(Int, pos)
            positions[][idx[]] = ipos
            notify(positions)
            return Consume(true)
        end
        return Consume(false)
    end
end

function pointnear(f, positions, ipos, accuracy)
    for i in eachindex(positions)[end:-1:1]
        p = positions[i]
        if p[1] in (ipos[1]-accuracy:ipos[1]+accuracy) && 
           p[2] in (ipos[2]-accuracy:ipos[2]+accuracy)
            # Remove a point
            return f(i)
            break
        end
    end
    return nothing
end

function ison(line, point, accuracy)
    (x1, y1), (x2, y2) = line
    x = round(point[1])
    y = round(point[2])
    grad = (y2 - y1) / (x2 - x1)
    if grad in (Inf, -Inf, NaN, NaN32)
        return x2 == x && inbounds((y1, y2), y)
    elseif grad == 0
        return y2 == y && inbounds((x1, x2), x)
    else
        inbounds((y1, y2), y) && inbounds((x1, x2), x) || return false
        if grad > -1 && grad < 1
            line_y = round(grad * (x - x1) + y1)
            return y in (line_y - accuracy):(line_y + accuracy)
        else
            line_x = round((y - y1)/grad + x1)
            return x in (line_x - accuracy):(line_x + accuracy)
        end
    end
end

inbounds((x1, x2), x) = x >= min(x1, x2) && x <= max(x1, x2)

function selectmultiple(A, fig, ax)
    if eltype(A) <: Colorant
        Makie.heatmap!(ax, Float64.(Gray.(A)))
    else
        Makie.heatmap!(ax, A)
    end
    positions = Observable(Point2{Float32}[])
    sct = Makie.scatter!(ax, positions, color=1:30, colormap=:reds)
    dragselect!(fig, ax, sct, positions, size(A))
    return positions
end

function _fitlinearmodels(points)
    x_model = lm(@formula(x_unknown ~ x_known), points)
    y_model = lm(@formula(y_unknown ~ y_known), points)
    return x_model, y_model
end

function linearwarp(A; template, points=nothing, models::Union{Nothing,Tuple}=nothing, missingval=missing)
    @show missingval
    x_model, y_model = if isnothing(models) 
        isnothing(points) && error("pass either `points::Tuple` to fit or fitted `models::Tuple`")
        _fitlinearmodels(points)
    else
        models
    end

    knownpoints = vec(collect((x_known = x, y_known=y) for (x, y) in Tuple.(CartesianIndices(template))))
    xs = round.(Int, predict(x_model, knownpoints))
    ys = round.(Int, predict(y_model, knownpoints))
    Awarped = Raster(similar(template, promote_type(typeof(missingval), eltype(A))); missingval)
    Awarped .= Rasters.missingval(Awarped)
    for (Ik, Iu) in  zip(CartesianIndices(template), CartesianIndex.(zip(xs, ys)))
        if checkbounds(Bool, A, Iu)
            Awarped[Ik] = A[Iu]
        end
    end
    return Awarped
end

end
