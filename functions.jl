using DynamicGrids
using DynamicGrids.Neighborhoods
using DynamicGrids.Neighborhoods: Window
using DataFrames

namedkeys(nt::NamedTuple{K}) where K = NamedTuple{K}(K)

abstract type SlopeFilter end
abstract type SlopeConvolution <: SlopeFilter end
# Not all slope algorithms can provide aspect
abstract type SlopeAspectConvolution <: SlopeConvolution end

# 
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


for f in (:slope, :aspect, :slopeaspect)
    f_filter = Symbol(f, :_filter) 
    @eval begin 
        function $(f)(elevation::Raster, method=FD2())
            padded = Neighborhoods.addpadding(parent(elevation), 1; padval=missingval(elevation))
            newdata = $(f)(padded, FD2())
            @show size(elevation) size(padded) size(newdata)
            rebuild(elevation; data=newdata, name=$(QuoteNode(f))) 
        end
        function $(f)(elevation::AbstractArray, method=FD2())
            window = Window{1}()
            Neighborhoods.broadcast_neighborhood(window, elevation) do w, e
                $(f_filter)(method, w, e)
            end
        end
    end
end

function nearest_distances(presences::AbstractArray{Bool}) 
    maxval = sqrt(sum(size(presences) .^ 2))
    nearest_distances!(fill(maxval, dims(presences)), presences)
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

function clean_categories(src::AbstractArray; categories=(), neighborhood=Moore{2,2}(), keep_neigborless=false, missingval=missing, despecle=true)
    counts = zeros(length(categories))
    ax = unpad_axes(src, neighborhood)
    dst = similar(src, promote_type(eltype(src), typeof(missingval)))
    dst .= missingval
    broadcast!(view(dst, ax...), CartesianIndices(ax)) do I
        DynamicGrids.Neighborhoods.applyneighborhood(neighborhood, src, I) do hood, v
            catcounts = map(categories) do c
                ds = distances(hood)
                acc = zero(1/first(ds))
                for i in 1:length(hood)
                    n = hood[i]
                    if n !== missingval && n === c 
                        acc += 1/ds[i]
                    end
                end
                return acc
            end
            if all(==(0), catcounts)
                if keep_neigborless
                    return v
                else
                    return missingval
                end
            end
            if despecle && !isequal(v, missingval) && v in categories
                if (count(==(v), skipmissing(hood)) > length(hood) / 10)
                    return v
                else
                    c = categories[findmax(catcounts)[2]]
                    return c
                end
            else
                return categories[findmax(catcounts)[2]]
            end
        end
    end
    return dst
end

function plot_lc(lc_raster::Raster)
    Plots.plot(lc_raster; 
        color=palette(:Paired_12, 12), 
        clims=(1, 13), size=(2000,2000),
        colorbar_ticks=map(Pair, 0:12, lc_categories),
    )
end

function plot_lc_makie(lc_raster::Raster)
    fig = Makie.Figure()
    ax, hm = Makie.heatmap(fig[1, 1], parent(parent(dims(lc_raster, X))), parent(parent(dims(lc_raster, Y))), parent(read(lc_raster)),
        colormap=cgrad(:cyclic_mygbm_30_95_c78_n256, 13, categorical=true), colorrange=(0, 13)
    )
    # ax.aspect = Makie.AxisAspect(1)
    # Makie.Colorbar(fig[1, 2], hm; 
        # ticks=(0:12, lc_categories),
    # )
    return fig
end

function rasterize_lc(template, shape_file, crs_file; res=nothing, categories)
    lc_shape = Shapefile.Table(shape_file)
    lc_crs = WellKnownText(readlines(crs_file)[1])
    lc_df = DataFrame(lc_shape)
    lc_raster = Raster(similar(template, Int32); missingval=typemin(Int32))
    lc_raster .= typemin(Int32)
    if !isnothing(res)
        lc_raster = read(resample(lc_raster, res; crs=lc_crs))
    end
    display(dims(lc_raster))
    # Order of rasterization matters?... (probably should calculate areas?)
    fillvals = [
        categories.No_Data,
        categories.Water,
        categories.Barren_land,
        categories.Forest,
        categories.Shrub_vegetation,
        categories.Herbaceaous_vegetation,
        categories.Mangrove,
        categories.Other_cropland,
        categories.Sugarcane,
        categories.Pasture,
        categories.Continuous_urban,
        categories.Discontinuous_urban,
    ]
    @show fillvals
    for fillval in fillvals
        rows = filter(x -> x.ocsol_num == fillval, lc_df)
        if length(rows.geometry) > 0
            fillname = first(eachrow(rows)).ocsol_name
            @show fillname fillval
            rasterize!(lc_raster, rows.geometry; fill=fillval)
        end
    end
    return lc_raster
end
