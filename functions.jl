using DynamicGrids
using DynamicGrids.Neighborhoods
using DynamicGrids.Neighborhoods: Window

abstract type SlopeFilter end
abstract type SlopeConvolution <: SlopeFilter end

# 
struct FD2 <: SlopeConvolution end
struct FD3Reciprocal <: SlopeConvolution end
struct FD3ReciprocalSquared <: SlopeConvolution end
struct FD3Linear <: SlopeConvolution end
struct FDFrame <: SlopeConvolution end
struct SimpleD <: SlopeConvolution end

@inline function slope_filter(method::SlopeConvolution, n::Window)
    fx, fy = _slope_conv(method, n)
    return atan(√(fx^2 + fy^2))
end

@inline function _slope_conv(::FD2, n::Window)
    fx = (n[6] - n[4]) / 2n[5]
    fy = (n[8] - n[2]) / 2n[5] 
    return fx, fy
end

@inline function _slope_conv(::FD3Reciprocal, n::Window)
    fx = (n[3] -n[1] + √(2(n[6] - n[4])) + n[9] - n[7]) / (4 + 2 * √(2)) * n[5]
    fy = (n[7] -n[1] + √(2(n[8] - n[2])) + n[9] - n[3]) / (4 + 2 * √(2)) * n[5]
    return fx, fy
end

@inline function _slope_conv(::FD3Linear, n::Window)
    fx = (n[3] - n[1] + n[6] - n[4] + n[9] - n[7]) / 6n[5]
    fy = (n[7] - n[1] + n[8] - n[2] + n[9] - n[3]) / 6n[5]
    return fx, fy
end

@inline function _slope_conv(::FD3ReciprocalSquared, n::Window)
    fx = (n[3] - n[1] + 2(n[6] - n[4]) + n[9] - n[7]) / 8n[5]
    fy = (n[7] - n[1] + 2(n[8] - n[2]) + n[9] - n[3]) / 8n[5]
    return fx, fy
end

@inline function _slope_conv(::FDFrame, n::Window)
    fx = (n[3] - n[1] + n[9] - n[7]) / 4n[5]
    fy = (n[7] - n[1] + n[9] - n[3]) / 4n[5] 
    return fx, fy
end

@inline function _slope_conv(::SimpleD, n::Window)
    fy = (n[5] - n[2]) / n[5]
    return fx, fy
end


struct MaxSlope <: SlopeFilter end

@inline function slope_filter(method::MaxSlope, n::Window)
    g = n[5] 
    slopes = (
        abs((n[5] - n[2]) / g), 
        abs((n[5] - n[4]) / g),
        abs((n[5] - n[6]) / g), 
        abs((n[5] - n[8]) / g), 
        abs((n[5] - n[1]) / (√(2)*g)), 
        abs((n[5] - n[3]) / (√(2)*g)), 
        abs((n[5] - n[7]) / (√(2)*g)), 
        abs((n[5] - n[9]) / (√(2)*g)), 
    )
    xmissing = map(ismissing, slopes)
    if all(xmissing)
        return missing
    elseif any(xmissing)
        return maximum(skipmissing(slopes))
    else
        return maximum(slopes)
    end
end

function slope(elevations::AbstractArray, method=FD2())
    window = Window{1}()
    Neighborhoods.broadcast_neighborhood(window, elevations) do w
        slope_filter(method, w)
    end
end

function nearest_distances(presences::AbstractArray{Bool}) 
    nearest_distances!(fill(Inf64, dims(presences)), presences)
end
function nearest_distances!(
    distances::AbstractArray{<:AbstractFloat}, presences::AbstractArray{Bool}
)
    distlookup = broadcast(CartesianIndices(size(distances))) do I
        sqrt(sum((Tuple(I) .- 1) .^ 2))
    end
    for I in CartesianIndices(presences) 
        # If the cell value is true, add distances to
        # it to the D array, otherwise skip
        presences[I] || continue
        for D in CartesianIndices(presences) 
            L = CartesianIndex(abs.(Tuple(D) .- Tuple(I)) .+ 1)
            if distances[D] > distlookup[L]
                 distances[D] = distlookup[L]  
            end
        end
    end
    return distances
end


function clean_categories(src::AbstractArray; categories=(), neighborhood=Window{2,2}(), missingval=missing)
    counts = zeros(length(categories))
    ax = unpad_axes(src, neighborhood)
    dst = similar(src, promote_type(eltype(src), typeof(missingval)))
    broadcast!(view(dst, ax...), view(src, ax...), CartesianIndices(ax)) do v, I
        DynamicGrids.Neighborhoods.applyneighborhood(neighborhood, src, I) do hood
            ismissing(v) && return v
            catcounts = map(categories) do c
                acc = zero(first(categories))
                for (n, d) in zip(neighbors(hood), distances(hood))
                    if !ismissing(n) && n == c 
                        acc += 1/d
                    end
                end
                return acc
            end
            if v in categories
                # de-speckle
                if v in hood 
                    categories[findmax(catcounts)[2]]
                else
                    @show "speckle found"
                    missingval
                end
            else
                if all(==(0), catcounts)
                    missingval
                else
                    categories[findmax(catcounts)[2]]
                end
            end
        end
    end
end

