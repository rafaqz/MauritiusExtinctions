using Rasters
using DimensionalData
using DynamicGrids
using DynamicGrids.Neighborhoods
using DynamicGrids.Neighborhoods: Window
using OffsetArrays

function _initialise_cost!(active, acc, origins, hood)
    for I in CartesianIndices(origins)
        if @inbounds origins[I] > zero(eltype(origins))
            @inbounds acc[I] = zero(eltype(acc))
            # Add an index to the `active` set if it has any neighbors not
            # in `origins` as we only need to use the edges of the origin areas.
            Neighborhoods.applyneighborhood(hood, origins, I) do hood, val
                if any(map(==(zero(eltype(origins))), hood))
                    push!(active, I)
                end
            end
        end
    end
    return active
end

"""
    cost_distance([f=meancost]; origins, costs, [missingval])

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
function cost_distance(f=meancost; costs, origins, kw...)
    # Create a grid we will update cost distance into.
    # The initial cost is the maximum `Float64` value. This will later
    # be the `missingval`, as it will remain in any cells that were not touched by
    # the algorithm (because one of the other arrays has missing values there).
    acc = OffsetArray(fill(typemax(Float64), size(origins) .+ 2), map(s -> 0:s+1, size(origins)))
    # Write the cost-distance to the accumulator grid
    cost_distance!(f, acc; origins, costs, kw...)
    # Remove the padding edge, and wrap the output as an AbstractRaster
    # if costs was one, updating the name and `missingval`.
    if costs isa AbstractRaster
        return rebuild(costs; data=Neighborhoods.unpad_view(acc, 1), name=:cost, missingval=typemax(Float64))
    else
        return Neighborhoods.unpad_view(acc, 1)
    end
end
function cost_distance!(f, acc;
    origins, costs, cellsize=1, missingval=missingval(costs),
)
    # The neighborood is a simple 3 * 3 moore neighborhood (ring)
    hood = Moore{1}()
    # The active cells are a `Set` of `CartesianIndices` linking 
    # to positions in acc/origins/costs.
    # "active" means assigned their lowest cost in the last round, 
    # but with neighbors whose cost has not been recalulated.
    active = Set{CartesianIndex{2}}()
    new_active = Set{CartesianIndex{2}}()
    # We add the first active cells from the `origins` array
    _initialise_cost!(active, acc, origins, hood)
    # And calculate how many active cells there are
    n_active = length(active)
    # Now loop while there are still active cells that are
    # reducing cost-distance somewhere on the grid.
    while n_active > 0
        for I in active 
            # Get the current accumulated cost
            @inbounds a = acc[I]
            # Get the cell cost (may be a NamedTuple if costs is a RasterStack)
            @inbounds cell_cost = costs[I]
            # Missing cells are skipped
            ismissingval(cell_cost, missingval) && continue
            # Loop over the neighborhood offsets and distances from center cell
            for (O, d) in zip(Neighborhoods.cartesian_offsets(hood), Neighborhoods.distances(hood))
                # Get the index of the neighboring cell
                NI = O + I
                # Out of bounds cells are skipped because our costs are not padded
                checkbounds(Bool, costs, NI) || continue
                @inbounds neighbor_cost = costs[NI]
                # Missing valued neighbors are skipped
                ismissingval(neighbor_cost, missingval) && continue
                @inbounds cur_cost = acc[NI]
                # Calculate the new cost by adding the cost to get to the
                # neighbor to the cost to get to the current cell
                new_cost = a + f(cell_cost, neighbor_cost, d * cellsize)
                # Update the costs if the new cost is lower than any previous path
                if new_cost < cur_cost
                    @inbounds acc[NI] = new_cost
                    push!(new_active, NI)
                end
            end
        end
        # Remove the active cells, their neighbors have been calculated and 
        # added to new_active, so they are done unless a cheaper path reaches 
        # them again later on.
        empty!(active)
        # Swap the active and new active Dicts for the next round of calculations
        active, new_active = new_active, active
        # Update how many cells we are calculating next round
        n_active = length(active)
    end
    return acc
end

ismissingval(val, missingval) = val === missingval
ismissingval(vals::NamedTuple, missingval) = any(map(val -> ismissingval(val, missingval), vals))
ismissingval(vals::NamedTuple, missingvals::NamedTuple) = any(map(ismissingval, vals, missingvals)) 

"""
    meancost(a, b, distance)

Calculate the cost of moving between cells with resistance `a` and `b` over distance `d`
"""
meancost(a, b, d) = d * (a + b) / 2

struct SlopeCost
    slopefactor::Float64
    distfactor::Float64
end
SlopeCost(; slopefactor=-3.5, distfactor=0.001) = SlopeCost(slopefactor, distfactor)

"""
    (sc::SlopeCost)(a, b, distance)

Calculate the cost of moving between elevation `a` to elevation `b` over distance `d`
"""
(sc::SlopeCost)(a, b, d) = 6 / exp(sc.slopefactor * (abs(a - b) / d + sc.distfactor)) * d
# (sc::SlopeCost)(a, b, d) = (abs(a - b) / d + sc.distfactor) * d

"""
    (sc::CombinedCosts)(a, b, distance)

Works with a RasterStack input and multiple functions
in a `NamedTuple` specifying which layers to use.
"""
struct CombinedCost{T<:NamedTuple,O}
    funcs::T
    op::O
end
function (cc::CombinedCost)(as, bs, d)
    mapreduce(cc.op, as, bs, cc.funcs) do a, b, f
        f(a, b, d)
    end
end

# slopecost(slope) = 6 / exp(-3.5 * (slope + 0.05))
