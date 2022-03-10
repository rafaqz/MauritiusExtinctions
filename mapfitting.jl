module RasterUtils

using Rasters
using GLM
using Makie
using GeometryBasics
using Colors
using Tables
using Observables
using ColorVectorSpace

using Makie
using Rasters.LookupArrays

function manualwarp(As...; template::Raster, points=nothing, missingval=missing)
    points = select_common_points(A1; template, points, missingval)
    applywarp(As...; template, points, missingval)
end

function applywarp(As...; template::Raster, points=nothing, missingval=missing)
    A1 = first(As)
    if A1 isa Raster
        As = map(A -> reorder(A, ForwardOrdered), As)
    end
    template = reorder(template, ForwardOrdered)
    warped = _warp_from_points(As, template, points, missingval)
    # Show updated heatmap
    # display(Makie.heatmap(map(parent, dims(first(warped)))..., parent(first(warped))))
    if length(warped) == 1
        display(Makie.image(parent(warped)))
        return first(warped)
    else
        return warped
    end
end

function _warp_from_points(As::Tuple, template, points, missingval)
    models = _fitlinearmodels(points)
    return map(A -> linearwarp(A; template, models, missingval), As)
end

select_common_points(A; template, kw...) = _select_common_points(A, template; kw...)

_select_common_points(A, template::Raster; kw...) = 
    _select_common_points(A, parent(reorder(template, ForwardOrdered)); kw...)
function _select_common_points(A, template::AbstractArray; points=nothing, missingval)
    # map(A -> size(A) == size(first(As)), As) || throw(ArgumentError("Intput raster sizes are not the same"))
    fig = Figure()
    ax1 = Makie.Axis(fig[1,1]; title="Source raster with known crs/resolution - `template` kw")
    ax2 = Makie.Axis(fig[1,2]; title="First raster with unknown crs/resolution")
    ax1.aspect = ax2.aspect = Makie.AxisAspect(1)
    dragging1 = Ref(false)
    dragging2 = Ref(false)
    knownpoints, unknownpoints = if !isnothing(points) && Tables.rowcount(points) > 0
        table2points(points)
    else
        Point2{Float32}[], Point2{Float32}[]
    end
    @show knownpoints unknownpoints
    knownpoints = selectmultiple(parent(template), fig, ax1; dragging=dragging1, points=knownpoints)
    unknownpoints = selectmultiple(A, fig, ax2; dragging=dragging2, points=unknownpoints)
    @show knownpoints unknownpoints
    finallimits = Ref{Any}(nothing)
    overlay = nothing
    lift(knownpoints, unknownpoints) do k, u
        (dragging1[] || dragging2[]) && return nothing # Dont update during drag
        len = min(length(k), length(u))
        (length(k) == length(u) && len >= 3) || return nothing
        points = points2table(k[1:len], u[1:len])
        warped = linearwarp(A; template, points, missingval)
        finallimits[] = ax1.finallimits
        if !isnothing(overlay)
            delete!(ax1, overlay)
        end
        overlay = _heatmap!(ax1, parent(warped); colormap=(:viridis, 0.2)) 
        ax1.finallimits = finallimits[]
        return nothing
    end
    screen = display(fig)
    println("Select points in rasters, then close the window")
    while screen.window_open[] 
        sleep(0.1)
    end
    length(knownpoints[]) == length(unknownpoints[]) || error("Number of selected points must be the same for each raster")
    return points2table(knownpoints[], unknownpoints[])
end

function _heatmap!(ax, A; colormap=:viridis, transparency=false) 
    if eltype(A) <: Colorant
        Makie.heatmap!(ax, Float64.(Gray.(A)); colormap, transparency)
    else
        Makie.heatmap!(ax, A; colormap, transparency)
    end
end

function points2table(knownpoints::Vector, unknownpoints::Vector)
    knowntable = (((x, y),) -> (x_known=Float64(x), y_known=Float64(y))).(knownpoints)
    unknowntable = (((x, y),) -> (x_unknown=Float64(x), y_unknown=Float64(y))).(unknownpoints)
    return merge.(knowntable, unknowntable)
end

function table2points(table)
    knownpoints = Point2{Float32}.(collect(zip(Tables.getcolumn(table, :x_known), Tables.getcolumn(table, :y_known))))
    unknownpoints = Point2{Float32}.(collect(zip(Tables.getcolumn(table, :x_unknown), Tables.getcolumn(table, :y_unknown))))
    return knownpoints, unknownpoints
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

function dragselect!(fig, ax, sct, positions, pixelsize; 
    selected=Ref(false), dragging=Ref(false), caninsert=false
)
    selected[] = false
    # Get pixel click accuracy from the size of visable heatmap.
    accuracy = lift(ax.finallimits) do fl
        round(Int, maximum(fl.widths) / 100)
    end
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
                                if online
                                    insert = true
                                    idx[] = i + 1
                                    insert!(positions[], i + 1, ipos)
                                    notify(positions)
                                    break
                                end
                                lastp = p
                            end
                        end
                        if !insert
                            push!(positions[], ipos)
                            idx[] = lastindex(positions[])
                            notify(positions)
                        end
                    end
                    dragging[] = true 
                    selected[] = true
                else
                    selected[] = false
                end
            elseif event.action == Mouse.release
                dragging[] = false
                notify(positions)
            end
        # Delete points with right click
        elseif event.button == Mouse.right
            if pos_px in ax.scene.px_area[]                    
                pointnear(positions[], ipos, accuracy[]) do i
                    isnothing(i) || deleteat!(positions[], i)
                    notify(positions)
                end
            end
            selected[] = false
        end
        @show selected
        return Consume(dragging[])
    end
    # Mouse drag event
    on(events(fig).mouseposition, priority = 2) do mp
        if dragging[]
            pos = Makie.mouseposition(ax.scene)
            ipos = round.(Int, pos)
            # Check for sync problems
            # if ipos in eachindex(positions[])
            positions[][idx[]] = ipos
            notify(positions)
            # end
            return Consume(true)
        end
        return Consume(false)
    end
    on(events(fig).keyboardbutton) do event
        if selected[] && event.action in (Keyboard.press, Keyboard.repeat)
            event.key == Keyboard.right  && _move(positions, idx[], (1, 0))
            event.key == Keyboard.up     && _move(positions, idx[], (0, 1))
            event.key == Keyboard.left   && _move(positions, idx[], (-1, 0))
            event.key == Keyboard.down   && _move(positions, idx[], (0, -1))
        end
        # Let the event reach other listeners
        return Consume(false)
    end
end

function _move(positions, i, dir)
    positions[][i] = positions[][i] .+ dir 
    notify(positions)
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

function selectmultiple(A, fig, ax; transparency=false, points, kw...)
    _heatmap!(ax, A; transparency) 
    positions = Observable(points)
    @show typeof(points)
    sct = Makie.scatter!(ax, positions, color=1:30, colormap=:reds)
    labels = lift(p -> string.(1:length(p)), positions)
    Makie.text!(ax, labels; position=positions)
    dragselect!(fig, ax, sct, positions, size(A); kw...)
    return positions
end

function _fitlinearmodels(points)
    x_model = lm(@formula(x_unknown ~ x_known + y_known), points)
    y_model = lm(@formula(y_unknown ~ y_known + x_known), points)
    return x_model, y_model
end

function linearwarp(A; template, points=nothing, models::Union{Nothing,Tuple}=nothing, missingval=missing)
    # @show points size(A) size(template)
    x_model, y_model = if isnothing(models) 
        isnothing(points) && error("pass either `points::Tuple` to fit or fitted `models::Tuple`")
        _fitlinearmodels(points)
    else
        models
    end
    # display(x_model); display(y_model)
    pixelpoints = vec(collect((x_known = x, y_known=y) for (x, y) in Tuple.(CartesianIndices(template))))
    xs = round.(Int, predict(x_model, pixelpoints))
    ys = round.(Int, predict(y_model, pixelpoints))
    T = promote_type(typeof(missingval), eltype(A))
    Awarped = similar(template, T)
    Awarped .= missingval
    for (Ik, Iu) in  zip(CartesianIndices(template), CartesianIndex.(zip(xs, ys)))
        if checkbounds(Bool, A, Iu)
            Awarped[Ik] = A[Iu]
        end
    end
    return Awarped
end

end
