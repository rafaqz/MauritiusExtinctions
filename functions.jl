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

function clean_categories(src::AbstractArray; categories=(), neighborhood=Moore{2,2}(), missingval=missing)
    counts = zeros(length(categories))
    ax = unpad_axes(src, neighborhood)
    dst = similar(src, promote_type(eltype(src), typeof(missingval)))
    dst .= missingval
    broadcast!(view(dst, ax...), view(src, ax...), CartesianIndices(ax)) do v, I
        DynamicGrids.Neighborhoods.applyneighborhood(neighborhood, src, I) do hood
            catcounts = map(categories) do c
                acc = 0
                for (n, d) in zip(neighbors(hood), distances(hood))
                    if !isequal(n, missingval) && n == c 
                        acc += 1/d
                    end
                end
                return acc
            end
            if all(==(0), catcounts)
                return missingval
            end
            if !isequal(v, missingval) && v in categories
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


lc_categories = [
  "No Data",
  "Continuous urban",
  "Disontinuous urban",
  "Forest",
  "Shrub vegetation",
  "Herbaceaous vegetation",
  "Mangrove",
  "Barren land",
  "Water",
  "Sugarcane",
  "Pasture",
  "",
  "Other cropland",
]

function plot_lc(lc_raster)
    plot(lc_raster; 
        color=palette(:Paired_12, 12), 
        clims=(1, 13), size=(2000,2000),
        colorbar_ticks=1:12,
    )
end

function plot_lc_makie(lc_raster)
    fig = Figure()
    ax, hm = Makie.heatmap(fig[1, 1], parent(dims(lc_raster, X)), parent(dims(lc_raster, Y)), parent(lc_raster),
        colormap=cgrad(:cyclic_mygbm_30_95_c78_n256, 13, categorical=true), colorrange=(0, 13)
    )
    ax.aspect = AxisAspect(1)
    Colorbar(fig[1, 2], hm; 
        ticks=(0:12, lc_categories),
    )
    return fig
end

function rasterize_lc(template, shape_file, crs_file; res=50)
    lc_shape = Shapefile.Table(shape_file)
    lc_crs = WellKnownText(readlines(crs_file)[1])
    lc_df = DataFrame(lc_shape)
    lc_raster = Raster(similar(template, Int32); missingval=typemin(Int32))
    lc_raster .= typemin(Int32)
    lc_raster = read(resample(lc_raster, res; crs=lc_crs))
    c = Dict(map(=>, lc_categories, 0:12))
    # Order of rasterization matters?... (probably should calculate areas?)
    fillvals = [
        c["No Data"],
        c["Water"],
        c["Barren land"],
        c["Forest"],
        c["Shrub vegetation"],
        c["Herbaceaous vegetation"],
        c["Mangrove"],
        c["Other cropland"],
        c["Sugarcane"],
        c["Pasture"],
        c["Continuous urban"],
        c["Disontinuous urban"],
    ]
    # @show fillvals
    for fillval in fillvals
        rows = filter(x -> x.ocsol_num == fillval, lc_df)
        if length(rows.geometry) > 0
            fillname = first(eachrow(rows)).ocsol_name
            @show fillval, fillname
            rasterize!(lc_raster, rows.geometry; fill=fillval)
        else
            @show fillval
        end
    end
    display(plot_lc_makie(lc_raster))
    @show typeof(parent(lc_raster))
    display(parent(lc_raster))
    return lc_raster
end
