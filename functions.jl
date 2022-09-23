
namedkeys(nt::NamedTuple{K}) where K = NamedTuple{K}(K)

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
        categories.Herbaceaous_vegetation,
        categories.Shrub_vegetation,
        categories.Barren_land,
        categories.Other_cropland,
        categories.Sugarcane,
        categories.Pasture,
        categories.Forest,
        categories.Mangrove,
        categories.Continuous_urban,
        categories.Discontinuous_urban,
    ]
    for fillval in fillvals
        rows = filter(x -> x.ocsol_num == fillval, lc_df)
        if length(rows.geometry) > 0
            fillname = first(eachrow(rows)).ocsol_name
            rasterize!(lc_raster, rows.geometry; fill=fillval)
        end
    end
    return lc_raster
end
