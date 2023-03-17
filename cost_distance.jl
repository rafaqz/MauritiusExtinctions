using Rasters
using DimensionalData
using Neighborhoods
using Neighborhoods: Window
using OffsetArrays
using Unitful

function _initialise_cost!(active, accumulated_costs, origins, hood)
    for I in CartesianIndices(origins)
        o = @inbounds origins[I]
        if !ismissing(o) && o < Inf * u"hr"
            @inbounds accumulated_costs[I] = o
            # Add an index to the `active` set if it has any neighbors not
            # in `origins` as we only need to use the edges of the origin areas.
            if any(x -> !ismissing(x) && x == Inf * u"hr", Neighborhoods.neighborhood(origins, I))
                push!(active, I)
            end
        end
    end
    return active
end

"""
    cost_distance([f=meancost]; origins, elevation, [missingval])

Calculate the cost-distance between cells in `origins` and
all other points in the `origins`/`costs`.

`origins` is a matrix of `Real` where values larger than `1` are used
as starting points with cost of zero.

`costs` is a matrix where values to pass to `a` and `b` in function `f`

The function `f` accept 3 parameters: `a`, `b`, and `distance`:
- `a`: The value in `costs` at the current location.
- `b`: The value in `costs` at the next location.
- `distance`: The distance between `a` and `b` locations.

By default, `f=meancost` so that `costs` taken as a grid of travel costs.
"""
function cost_distance(origins, elevation, resistance=nothing; kw...)
    # Create a grid we will update cost distance into.
    # The initial cost is the maximum `Float64` value. This will later
    # be the `missingval`, as it will remain in any cells that were not touched by
    # the algorithm (because one of the other arrays has missing values there).
    hood = Moore{1}()
    accumulated_time = NeighborhoodArray(fill(Inf * u"hr", size(origins)), hood; padding=Conditional(), boundary=Remove(typemax(Float64)))
    origins = NeighborhoodArray(parent(origins), hood; padding=Conditional(), boundary=Remove(typemax(Float64)))
    # Write the cost-distance to the accumulator grid
    cost_distance!(accumulated_time, origins, elevation, resistance; hood, kw...)
    # Remove the padding edge, and wrap the output as an AbstractRaster
    # if costs was one, updating the name and `missingval`.
    if elevation isa AbstractRaster
        return rebuild(elevation; data=Neighborhoods.parent(accumulated_time), name=:time, missingval=Inf * u"hr")
    else
        return Neighborhoods.unpad_view(acc, 1)
    end
end
function cost_distance!(accumulated_time, origins, elevation, resistance=nothing; 
    cellsize=1, missingval=missingval(elevation), hood
)
    # The neighborood is a simple 3 * 3 moore neighborhood (ring)
    # The active cells are a `Set` of `CartesianIndices` linking
    # to positions in acc/origins/costs.
    # "active" means assigned their lowest cost in the last round,
    # but with neighbors whose cost has not been recalulated.
    active_cells = Set{CartesianIndex{2}}()
    new_active_cells = Set{CartesianIndex{2}}()
    # We add the first active cells from the `origins` array
    _initialise_cost!(active_cells, accumulated_time, origins, hood)
    # And calculate how many active cells there are
    n_active_cells = length(active_cells)
    # Now loop while there are still active cells that are
    # reducing cost-distance somewhere on the grid.
    while n_active_cells > 0
        for I in active_cells
            # Get the current accumulated cost
            @inbounds current_time = accumulated_time[I]
            # Missing cells are skipped
            e1 = elevation[I]
            ismissingval(e1, missingval) && continue
            # Loop over the neighborhood offsets and distances from center cell
            for (O, d) in zip(Neighborhoods.cartesian_offsets(hood), Neighborhoods.distances(hood))
                # Get the index of the neighboring cell
                NI = O + I
                checkbounds(Bool, elevation, NI) || continue
                e2 = elevation[NI]
                # Out of bounds cells are skipped because our costs are not padded
                # Missing valued neighbors are skipped
                ismissing(e2) && continue
                @inbounds current_neighbor_time = accumulated_time[NI]
                # Calculate the new time by adding the time to get to the
                # neighbor to the time to get to the current cell
                r = isnothing(resistance) ? 1 : resistance[I]
                comb = time_taken(e1, e2, d * cellsize, r)
                # @show comb
                new_neighbor_time = current_time + comb
                # @show new_neighbor_time current_neighbor_time
                # Update the time if the new time is lower than any previous path
                if new_neighbor_time < current_neighbor_time
                    @inbounds accumulated_time[NI] = new_neighbor_time
                    push!(new_active_cells, NI)
                end
            end
        end
        # Remove the active cells, their neighbors have been calculated and
        # added to new_active_cells, so they are done unless a cheaper path reaches
        # them again later on.
        empty!(active_cells)
        # Swap the active and new active Dicts for the next round of calculations
        active_cells, new_active_cells = new_active_cells, active_cells
        # Update how many cells we are calculating next round
        n_active_cells = length(active_cells)
    end
    return accumulated_time
end

ismissingval(val, missingval) = val === missingval
ismissingval(vals::NamedTuple, missingval) = any(map(val -> ismissingval(val, missingval), vals))
ismissingval(vals::NamedTuple, missingvals::NamedTuple) = any(map(ismissingval, vals, missingvals))

"""
    meancost(a, b, distance)

Calculate the cost of moving between cells with resistance `a` and `b` over distance `d`
"""
meancost(a, b, d) = d * (a + b) / 2

"""
    time_taken(e1, e2, d, v)

Calculate the cost of moving between elevation `e1` to elevation `e2`
over distance `d` through resistance `r`

cost = S * I * A * E
"""
function time_taken(e1, e2, d, r)
    s = walking_speed(ImhofTobler(), abs(e1 - e2), d)
    return d / (s * r)
end
# (sc::SlopeCost)(a, b, d) = (abs(a - b) / d + sc.distfactor) * d

abstract type AnisotropicCost end

struct ImhofTobler <: AnisotropicCost
    slopefactor::Float64
    distfactor::Float64
end
ImhofTobler(; slopefactor=-3.5, distfactor=0.05) =
    ImhofTobler(slopefactor, distfactor)

walking_speed(x::ImhofTobler, dvert, dhoriz) =
    walking_speed(x, dvert / dhoriz)
function walking_speed(x::ImhofTobler, slope)
    (6ℯ^(x.slopefactor * (slope + x.distfactor))) * u"km/hr"
end

# struct JabeWhite <: AnisotropicCost end

# anisotropic_cost(x::JabeWhite, dh, dx) = anisotropic_cost(x, dh / dx)
# function anisotropic_cost(::JabeWhite, S)
#     if S < -60
#         Inf
#     elseif S < -6
#         20.9tan(S)^2 + 4.18tan(S) + 1.38
#     elseif S < 60
#         52.1 * 10^3tan(S)^2 + 10.4tan(S) + 2.65
#     else
#         Inf
#     end
# end

abstract type IsotropicCost end

struct SouleGoldman <: IsotropicCost end

function land_cover_speed(::SouleGoldman, category)
    relative_movement_speed[category]
end

const movement_speed = (
    road=5.0u"km/h",
    dirt_road=4.0u"km/h",
    track=3.2u"km/h",
    cleared_land=2.0u"km/h",
    forest=1.0u"km/h",
    dense_forest=0.5u"km/h",
    wetland=0.5u"km/h",
)

const relative_movement_speed = map(movement_speed) do s
    s / movement_speed.road
end

