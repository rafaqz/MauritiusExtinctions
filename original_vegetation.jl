using Rasters, GLMakie, ColorSchemes, Extents
using DBFTables
# using CairoMakie
GLMakie.activate!()
include("landcover_compilation.jl")

function plot_habitats!(fig, data; 
    colormap, nrows, ncols, show_uncertain=true
) 
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
            titlesize=20,
        )
        tight_ticklabel_spacing!(ax)
        if show_uncertain
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
            poly!(ax, rect; color=stripe, strokewidth=0)
            Makie.heatmap!(ax, stripemask; colormap=whites, colorrange=(0, 1))
        end
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
        lines!(ax, l, a; color, ticksize=14)
        band!(ax, l, base, a; color)
        pa = Point2f.(l, a)
        pb = Point2f.(l, b)
        polygon = [pa..., pb[end:-1:1]..., pa[1]]
        poly!(ax, polygon; color=stripe, strokewidth=0)#, alpha=0.5)
        base = b
    end
end

function _legend!(position, habitat_colors, habitat_names)
    fig = position.layout.parent 
    # Legend
    habitat_elements = map(habitat_colors) do color
        PolyElement(; color, strokewidth=0)
    end
    names = replace.(habitat_names, Ref('_' => ' '))
    Legend(position, habitat_elements, names, "Habitat class"; 
        titlesize=22,
        framevisible=false,
        labelsize=16,
        patchsize=(30.0f0, 30.0f0)
    )
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
fig = Figure(; size=(1700, 2000));
data = uncleared.mus
nrows = 5
ncols = 4
cmap = :tableau_20
habitat_colors = map(x -> getproperty(ColorSchemes, cmap)[(x - 1) / 9 ], 1:nhabitats.mus) |> reverse
# Heatmaps
plot_habitats!(fig, data; colormap=habitat_colors, nrows, ncols, show_uncertain=true) 
# Legend
_legend!(fig[5, 4], habitat_colors, island_habitat_names.mus)
# Area plot
line_ax = Axis(fig[6, 1:4])
plot_aggregate!(line_ax, data, habitat_colors)
# Pad the line plot little
rowgap!(fig.layout, 5, 40)
# Title
fig[7, :] = Label(fig, "Mauritius Habitat Loss (striped/transparent areas uncertain)"; fontsize=30)
save("images/mauritius_habitat_loss.png", fig)
# display(fig)

# Reunion
fig = Figure(; size=(1700, 2000));
data = uncleared.reu
nrows = 4
ncols = 3
cmap = :tableau_20
habitat_colors = map(x -> getproperty(ColorSchemes, cmap)[(x - 1) / (nhabitats.reu - 1) ], 1:nhabitats.reu)
# Heatmaps
plot_habitats!(fig, data; colormap=habitat_colors, nrows, ncols, show_uncertain=true);
_legend!(fig[3:4, 3], habitat_colors, island_habitat_names.reu)
# Area plot
line_ax = Axis(fig[5, 1:3])
plot_aggregate!(line_ax, data, habitat_colors)
# Pad the line plot little
rowgap!(fig.layout, 4, 100)
# Title
fig[6, :] = Label(fig, "Reunion Habitat Loss (striped/transparent areas uncertain)"; fontsize=30)
save("images/reunion_habitat_loss.png", fig)
# display(fig)

p = Makie.heatmap(uncleared.mus.certain[Ti=End-2])
Makie.heatmap!(p.axis, uncleared.mus.uncertain[Ti=End-2]; alpha=0.3)
Makie.heatmap!(p.axis, native_veg.mus; alpha=0.7, colormap=:reds)

p = Makie.heatmap(uncleared.reu.certain[Ti=End])
Makie.heatmap!(p.axis, uncleared.reu.uncertain[Ti=End]; alph=0.5)
Makie.heatmap!(p.axis, native_veg.reu; alpha=0.5, colormap=:reds)
