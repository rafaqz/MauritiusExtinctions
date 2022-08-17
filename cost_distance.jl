using DimensionalData
using DynamicGrids
using DynamicGrids.Neighborhoods
using DynamicGrids.Neighborhoods: Window
using OffsetArrays

function _initialise_cost!(active, acc, origins, hood)
    for I in CartesianIndices(origins)
        if @inbounds origins[I] > zero(eltype(origins))
            # Add index if it has any neighbors not in `origins`
            Neighborhoods.apply_neighborhood(hood, origins, I) do hood, val
                @inbounds acc[I] = zero(eltype(acc))
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
    acc = OffsetArray(fill(typemax(Float64), size(origins) .+ 2), map(s -> 0:s+1, size(origins)))
    cost_distance!(f, acc; origins, costs, kw...)
end
function cost_distance!(f, acc;
    origins, costs, cellsize=1, checkmissing=isequal(missingval(costs)),
)
    hood = Moore{1}()
    active = Set{CartesianIndex{2}}()
    new_active = Set{CartesianIndex{2}}()
    _initialise_cost!(active, acc, origins, hood)
    n_active = length(active)
    while n_active > 0
        for I in active 
            @inbounds a = acc[I]
            @inbounds cell_cost = costs[I]
            checkmissing(cell_cost) && continue
            for (O, d) in zip(Neighborhoods.cartesian_offsets(hood), Neighborhoods.distances(hood))
                NI = O + I
                checkbounds(Bool, costs, NI) || continue
                @inbounds neighbor_cost = costs[NI]
                checkmissing(neighbor_cost) && continue
                @inbounds cur_cost = acc[NI]
                new_cost = a + f(cell_cost, neighbor_cost, d * cellsize)
                if cur_cost > new_cost
                    @inbounds acc[NI] = new_cost
                    push!(new_active, NI)
                end
            end
        end
        empty!(active)
        active, new_active = new_active, active
        n_active = length(active)
    end
    if costs isa AbstractDimArray
        return rebuild(costs; data=Neighborhoods.unpad_view(acc, 1), name=:cost, missingval=typemax(Float64))
    else
        return Neighborhoods.unpad_view(acc, 1)
    end
end

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
