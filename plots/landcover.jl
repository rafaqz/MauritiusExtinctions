using GLMakie
using ColorSchemes

include("landcover_compilation.jl")
batlow = map(1:6) do i
    ColorSchemes.batlow[(i - 1) / 5]
end

function _nt_vecs(xs::AbstractArray{<:NamedVector{K}}) where K
    xs = dropdims(xs; dims=(X(), Y()))
    vecs = map(1:length(K)) do i
        getindex.(xs, i)
    end
    RasterStack(vecs...; name=K)
end

function plot_timeline(timeline_counts, striped, npixels; show=keys(timeline_counts))
    l = lookup(first(timeline_counts), Ti)
    xticks = eachindex(l)
    xtickformat = i -> string.(getindex.(Ref(l), Int.(i)))
    x = eachindex(l)
    fig = Figure(size=(2000, 2000));#, backgroundcolor="#a5b4b5")
    all_axes = map(enumerate(show)) do (j, statistic)
        heatmap_axes = map(1:length(x)) do i 
            A = striped[statistic][Ti=i]
            ax = Axis(fig[j*2 - 1, i]; autolimitaspect=1)
            tight_ticklabel_spacing!(ax)
            Makie.image!(ax, A; colormap=:batlow, colorrange=(1, 6), interpolate=false)
            hidedecorations!(ax)
            hidespines!(ax)
            ax
        end
        line_axis = Axis(fig[j*2, 1:length(x)];
            backgroundcolor=:white, 
            ylabel=titlecase(string(statistic)), 
            limits=((first(xticks) - 0.5, last(xticks) + 0.5), nothing),
            xticks, xtickformat,
        )
        line_axis.xzoomlock = true
        if statistic != :merged 
            hidexdecorations!(line_axis)
        end
        hidespines!(line_axis)
        for (i, k) in enumerate(classnames)
            if statistic == :merged
                y = timeline_counts[:merged][k] ./ npixels
                z = (timeline_counts[:merged][k] .- timeline_counts[:uncertain][k]) ./ npixels
                Makie.band!(line_axis, x, y, z; color=batlow[i], alpha=0.4)
                Makie.lines!(line_axis, x, z; color=batlow[i], linewidth=2)
            else
                y = timeline_counts[statistic][k] ./ npixels
                Makie.lines!(line_axis, x, y; color=batlow[i], linewidth=2)
            end
        end
        heatmap_axes, line_axis
    end
    Makie.linkaxes!(Iterators.map(last, all_axes)...)
    Makie.linkaxes!(Iterators.flatten(Iterators.map(first, all_axes))...)
    colgap!(fig.layout, Relative(0.001))
    rowgap!(fig.layout, Relative(0.001))
    return fig
end

timeline_counts = map(statistics) do s
    (c, f, u, a, re, m) = s
    combined = sum(+, c; dims=(X, Y)) |> _nt_vecs
    filled = sum(+, f; dims=(X, Y)) |> _nt_vecs
    uncertain = sum(+, u; dims=(X, Y)) |> _nt_vecs
    added = sum(+, a; dims=(X, Y)) |> _nt_vecs
    removed = sum(+, re; dims=(X, Y)) |> _nt_vecs
    merged = sum(+, m; dims=(X, Y)) |> _nt_vecs
    (; combined, filled, uncertain, added, removed, merged)
end
striped_statistics = stripe_raster(statistics, states)
classnames = keys(timeline_counts.mus.combined)
npixels = map(count, masks)

fig = plot_timeline(timeline_counts.rod, striped_statistics.rod, npixels.rod)
fig = plot_timeline(timeline_counts.reu, striped_statistics.reu, npixels.reu)
fig = plot_timeline(timeline_counts.mus, striped_statistics.mus, npixels.mus)
fig = plot_timeline(timeline_counts.rod, striped_statistics.rod, npixels.rod; show=(:merged,))
fig = plot_timeline(timeline_counts.reu, striped_statistics.reu, npixels.reu; show=(:merged,))
fig = plot_timeline(timeline_counts.mus, striped_statistics.mus, npixels.mus; show=(:merged,))

# save("images/mus_map_timeline.png", fig)
# display(fig)


