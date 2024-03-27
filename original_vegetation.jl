using Rasters, GLMakie, ColorSchemes, Extents
using DBFTables
using CairoMakie
includet("raster_common.jl")
mus_veg_path = "/home/raf/PhD/Mascarenes/Data/Selected/Mauritius/Undigitised/page33_mauritius_vegetation_colored.tif"
reu_veg_path = "/home/raf/PhD/Mascarenes/Data/Dominique/Vegetation_Rasters/pastveg3.tif"

# Build auxiliary rasters
# lc_predictions_paths = (
#     mus="$outputdir/lc_predictions_mus.nc",
#     reu="$outputdir/lc_predictions_reu.nc",
#     rod="$outputdir/lc_predictions_rod.nc",
# )
# original_veg = (;
#     mus=reorder(replace_missing(Raster(mus_veg_path), 0), lc_predictions.mus),
#     reu=reorder(resample(replace_missing(Raster(reu_veg_path), 0); to=masks.reu), lc_predictions.reu),
# )
# # netcdf has the annoying center locus for time
# lc_predictions = map(lc_predictions_paths) do path
#     lc_predictions = RasterStack(path) |>
#         x -> rebuild(Rasters.modify(BitArray, x); missingval=false) |>
#         x -> Rasters.set(x, Ti => Int.(maybeshiftlocus(Start(), dims(x, Ti), )))
# end
# original_veg.mus
# lc_predictions.mus

# island_veg_change = map(lc_predictions[keys(original_veg)], original_veg) do p, v
#     rebuild(UInt8.(broadcast_dims(*, p.native, v)); missingval=0)
# end
# p = Rasters.rplot(island_veg_change.mus[Ti=At(1700:2018)]; colormap=:viridis)
# save("images/mus_original_veg.png", p)


# habitat_sums = map(island_veg_change, island_habitat_names, nhabitats) do veg_change, habitat_names, nhabitat
#     As = map(1:nhabitat) do habitat
#         dropdims(sum(==(habitat), veg_change; dims=(X, Y)); dims=(X, Y))
#     end
#     cat(As...; dims=Dim{:habitat}(habitat_names))
# end
# k = :reu
# k = :mus
# cum = cumsum(habitat_sums[k]; dims=2)
# x = lookup(habitat_sums[k], Ti)
# fig = Figure()
# ax = Axis(fig[1, 1])
# for i in nhabitats[k]:-1:1
#     y = parent(cum[habitat=i])
#     Makie.lines!(x, y; color=:black)
#     band!(x, fill(0, length(x)), y; label = "Label")
# end
# fig[1, 2] = Legend(fig, ax, habitat_names)

function plot_habitats!(fig, data; colormap, nrows, ncols) 
    whites = [RGB(1), RGB(1)] 
    axs = map(axes(data.certain, Ti)) do i
        stripe = Makie.LinePattern(; 
            direction=Vec2f(1), width=5, tilesize=(20, 20),
            linecolor=(:grey, 0.7), background_color=(:white, 0.0)
        )
        n = length(axes(data.certain, Ti))
        r = rem(n, i)
        ax = Axis(fig[reverse(fldmod1(i, nrows))...]; 
            # aspect=DataAspect(),
            autolimitaspect=1,
            title=string(lookup(data.certain, Ti)[i]),
        )
        tight_ticklabel_spacing!(ax)
        uncertain = data.uncertain[Ti=i]
        Makie.heatmap!(ax, uncertain; alpha=0.5, colormap)
        stripemask = map(uncertain) do x
            x > 0 ? missing : 1
        end
        bs = Rasters.bounds(stripemask)
        rect = Polygon([
            Point2f(bs[1][1], bs[2][1]), 
            Point2f(bs[1][1], bs[2][2]), 
            Point2f(bs[1][2], bs[2][2]), 
            Point2f(bs[1][2], bs[2][1]), 
            Point2f(bs[1][1], bs[2][1]), 
        ])
        Makie.poly(rect)
        poly!(ax, rect; color=stripe, strokewidth=0)
        Makie.heatmap!(ax, stripemask; colormap=whites, colorrange=(0, 1))
        Makie.heatmap!(ax, data.certain[Ti=i]; colormap)
        hidedecorations!(ax)
        hidespines!(ax)
        ax
    end
    linkaxes!(axs...)
    axs
end

function plot_aggregate!(ax, data, habitat_colors)
    npixels = count(>(0), view(data.certain, Ti=1))
    certain_agg = map(eachindex(habitat_colors)) do i 
        dropdims(count(data.certain; dims=(X, Y)) do x
            x == i
        end; dims=(X, Y))
    end
    uncertain_agg = map(eachindex(habitat_colors)) do i 
        dropdims(count(data.uncertain; dims=(X, Y)) do x
            x == i
        end; dims=(X, Y))
    end
    # hidedecorations!(line_ax)
    hidespines!(ax)
    base = map(_ -> 0.0, certain_agg[1])

    for i in reverse(eachindex(habitat_colors))
        color = habitat_colors[i]
        stripe = Makie.LinePattern(; 
            direction=Vec2f(1), width=5, tilesize=(20, 20),
            linecolor=(:grey, 0.7), background_color=(color, 0.7)
        )
        a = certain_agg[i] ./ npixels .+ base
        b = (certain_agg[i] .+ uncertain_agg[i]) ./ npixels .+ base
        l = parent(lookup(a, Ti))
        lines!(ax, l, a; color)
        band!(ax, l, base, a; color)
        pa = Point2f.(l, a)
        pb = Point2f.(l, b)
        polygon = [pa..., pb[end:-1:1]..., pa[1]]
        poly!(ax, polygon; color=stripe, strokewidth=0)#, alpha=0.5)
        base = b
    end
end

certain_uncleared = map(statistics) do island
    map(island.merged) do xs
        count(xs) == 1 && xs.native
    end
end
uncertain_uncleared = map(statistics) do island
    map(island.merged) do xs
        count(xs) > 1 && xs.native
    end
end
original_veg = (;
    mus=reorder(replace_missing(Raster(mus_veg_path), 0), uncertain_uncleared.mus),
    reu=reorder(resample(replace_missing(Raster(reu_veg_path), 0); to=masks.reu), uncertain_uncleared.reu),
)
k = keys(original_veg)
uncleared = map(original_veg, certain_uncleared[k], uncertain_uncleared[k]) do v, c, u
    (; certain=broadcast_dims(*, c, v), uncertain=broadcast_dims(*, u, v))
end
reu_veg_df = DBFTables.Table("/home/raf/PhD/Mascarenes/Data/Dominique/Vegetation_Rasters/pastveg3.tif.vat.dbf") |> DataFrame
island_habitat_names = (;
    mus=["semi-dry_evergreen_forest", "open_dry_palm-rich_woodland", "wet_forest", "pandanus_swamp", "mossy_rainforest", "mangrove", "wetland vegetation"],
    reu=reu_veg_df.HABITAT,
)
nhabitats = map(length, island_habitat_names)

# Mauritius
fig = Figure(; size=(2000, 2400));
data = uncleared.mus
nrows = 5
ncols = 4
cmap = :tableau_20
habitat_colors = map(x -> getproperty(ColorSchemes, cmap)[(x - 1) / 9 ], 1:nhabitats.mus) |> reverse
# Heatmaps
plot_habitats!(fig, data; colormap=habitat_colors, nrows, ncols) 
# Legend
habitat_elements = map(habitat_colors) do color
    PolyElement(; color, strokewidth=0)
end
habitat_names = titlecase.(replace.(island_habitat_names.mus, Ref('_' => ' ')))
fig[5, 4] = Legend(fig, habitat_elements, habitat_names, "Habitat classes"; framevisible=false)
# Area plot
line_ax = Axis(fig[6, 1:4])
plot_aggregate!(line_ax, data, habitat_colors)
# Pad the line plot little
rowgap!(fig.layout, 5, 40)
# Title
fig[7, :] = Label(fig, "Mauritius Habitat Loss (striped/transparent areas uncertain)"; fontsize=30)
display(fig)
save("images/mauritius_habitat_loss.png", fig)

# Reunion
fig = Figure(; size=(2000, 2400));
data = uncleared.reu
nrows = 4
ncols = 3
cmap = :tableau_20
habitat_colors = map(x -> getproperty(ColorSchemes, cmap)[(x - 1) / (nhabitats.reu - 1) ], 1:nhabitats.reu)
# Heatmaps
plot_habitats!(fig, data; colormap=habitat_colors, nrows, ncols);
# Legend
habitat_elements = map(habitat_colors) do color
    PolyElement(; color, strokewidth=0)
end
habitat_names = replace.(island_habitat_names.reu, Ref('_' => ' '))
fig[3:4, 3] = Legend(fig, habitat_elements, habitat_names, "Habitat classes"; 
    framevisible=false,
    labelsize=20,
    patchsize=(30.0f0, 30.0f0)
)
# Area plot
line_ax = Axis(fig[5, 1:3])
plot_aggregate!(line_ax, data, habitat_colors)
# Pad the line plot little
rowgap!(fig.layout, 4, 40)
# Title
fig[6, :] = Label(fig, "Reunion Habitat Loss (striped/transparent areas uncertain)"; fontsize=30)
display(fig)
save("images/reunion_habitat_loss.png", fig)
