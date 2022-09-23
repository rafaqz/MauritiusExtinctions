using DynamicGrids
using DynamicGrids.Neighborhoods
using DynamicGrids.Neighborhoods: Window
using DataFrames

abstract type SlopeFilter end
abstract type SlopeConvolution <: SlopeFilter end
# Not all slope algorithms can provide aspect
abstract type SlopeAspectConvolution <: SlopeConvolution end

struct FD2 <: SlopeAspectConvolution end
struct FD3Reciprocal <: SlopeAspectConvolution end
struct FD3ReciprocalSquared <: SlopeAspectConvolution end
struct FD3Linear <: SlopeAspectConvolution end
struct FDFrame <: SlopeAspectConvolution end
struct SimpleDifference <: SlopeConvolution end

@inline function aspect_filter(method::SlopeConvolution, n::Window, val)
    fx, fy = _slope_conv(method, n, val)
    return _aspect(fx, fy)
end

@inline function slope_filter(method::SlopeConvolution, n::Window, val)
    fx, fy = _slope_conv(method, n, val)
    return _slope(fx, fy)
end

@inline function slopeaspect_filter(method::SlopeAspectConvolution, n::Window, val)
    fx, fy = _slope_conv(method, n, val)
    return _slope(fx, fy), _aspect(fx, fy)
end

_slope(fx, fy) = atan(√(fx^2 + fy^2))
function _aspect(fx, fy)
    (ismissing(fx) || ismissing(fy)) && return missing
    # Rotate - we want high Y (north) as the origin
    # TODO: pass through the Order for X/Y dims 
    # So the result always has zero at North
    -atan(fx, fy) 
end

@inline function _slope_conv(::FD2, n::Window, val)
    fx = (n[6] - n[4]) / 2val
    fy = (n[8] - n[2]) / 2val
    return fx, fy
end

@inline function _slope_conv(::FD3Reciprocal, n::Window, val)
    fx = (n[3] -n[1] + √(2(n[6] - n[4])) + n[9] - n[7]) / (4 + 2 * √(2)) * val
    fy = (n[7] -n[1] + √(2(n[8] - n[2])) + n[9] - n[3]) / (4 + 2 * √(2)) * val
    return fx, fy
end

@inline function _slope_conv(::FD3Linear, n::Window, val)
    fx = (n[3] - n[1] + n[6] - n[4] + n[9] - n[7]) / 6val
    fy = (n[7] - n[1] + n[8] - n[2] + n[9] - n[3]) / 6val
    return fx, fy
end

@inline function _slope_conv(::FD3ReciprocalSquared, n::Window, val)
    fx = (n[3] - n[1] + 2(n[6] - n[4]) + n[9] - n[7]) / 8val
    fy = (n[7] - n[1] + 2(n[8] - n[2]) + n[9] - n[3]) / 8val
    return fx, fy
end

@inline function _slope_conv(::FDFrame, n::Window, val)
    fx = (n[3] - n[1] + n[9] - n[7]) / 4val
    fy = (n[7] - n[1] + n[9] - n[3]) / 4val
    return fx, fy
end

@inline function _slope_conv(::SimpleDifference, n::Window, val)
    fy = (n[5] - n[2]) / n[5]
    return fx, fy
end


struct MaxSlope <: SlopeFilter end

@inline function slope_filter(method::MaxSlope, n::Window, g)
    # slopes = (
    #     abs((g - n[2]) / g), 
    #     abs((g - n[4]) / g),
    #     abs((g - n[6]) / g), 
    #     abs((g - n[8]) / g), 
    #     abs((g - n[1]) / (√(2)*g)), 
    #     abs((g - n[3]) / (√(2)*g)), 
    #     abs((g - n[7]) / (√(2)*g)), 
    #     abs((g - n[9]) / (√(2)*g)), 
    # )
    # xmissing = map(ismissing, slopes)
    # if all(xmissing)
    #     return missing
    # elseif any(xmissing)
    #     return maximum(skipmissing(slopes))
    # else
    #     return maximum(slopes)
    # end
end


for (f, filt) in (:slope => :slope_filter, :aspect => :aspect_filter, :_slopeaspect => :slopeaspect_filter)
    @eval begin 
        function $(f)(elevation::Raster, method=FD2())
            padded = Neighborhoods.pad_array(parent(elevation), 1; padval=missingval(elevation))
            newdata = $(f)(padded, FD2())
            rebuild(elevation; data=newdata, name=$(QuoteNode(f))) 
        end
        function $(f)(elevation::AbstractArray, method=FD2())
            window = Window{1}()
            Neighborhoods.broadcast_neighborhood(window, elevation) do w, e
                $(filt)(method, w, e)
            end
        end
    end
end

function slopeaspect(elevation, method=FD2())
    sa = _slopeaspect(elevation, method)
    slope = first.(sa)
    aspect = last.(sa)
    nt = (; slope, aspect)
    return elevation isa Raster ? RasterStack(nt) : nt
end

