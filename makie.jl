using ColorSchemes

Makie.set_theme!(theme_light())

const COLORMAPS = [:magma, :viridis, :cividis, :inferno, :delta, :seaborn_icefire_gradient, :seaborn_rocket_gradient, :hot]

function mk(init, ruleset; maxpops=zero(eltype(init.pred_pop)), landcover, tspan, kw...)
    MakieOutput(init;
        kw...,
        tspan,
        fps=100,
        store=false,
        ruleset,
        sim_kw=(; printframe=true),
    ) do (; layout, frame, time)

        colorrange_obs = map(x -> Observable((zero(x), oneunit(x))), maxpops)
        # Landcover
        lc = lift(time) do i
            replace_missing(landcover[Ti(Near(tspan[i]))], NaN32)
        end
        ax_lc = Axis(layout[1, 1]; title="Landcover")
        # hidexdecorations!(ax_lc; grid=false)
        # hideydecorations!(ax_lc; grid=false)
        Makie.image!(ax_lc, lc; colormap=:batlow, colorrange=(0, 6), interpolate=false)

        # Predators
        pred_keys = propertynames(frame[].pred_pop[1])
        npreds = length(pred_keys)
        ncols = npreds + 1
        pred_axes = map(2:ncols) do i
            Axis(layout[1, i]; title=_title(pred_keys[i-1]))
        end
        hidexdecorations!.(pred_axes; grid=false)
        hideydecorations!.(pred_axes; grid=false)
        predators = map(1:npreds) do i
            Observable(rebuild(init.pred_pop, (x -> iszero(x) ? NaN : Float64(x)).(getindex.(frame[].pred_pop, i))))
        end
        on(frame) do f
            foreach(maximum(f.pred_pop), colorrange_obs) do m, obs 
                obs[] = (obs[][1], max(m, obs[][2]))
                notify(obs)
            end
            foreach(predators, 1:npreds) do pred, i
                pred[] .= (x -> iszero(x) ? NaN : Float64(x)).(getindex.(frame[].pred_pop, i))
                notify(pred)
            end
        end
        foreach(2:ncols, pred_axes, predators, pred_keys, COLORMAPS[1:npreds], colorrange_obs) do i, ax, pred, k, colormap, cr
            p = Makie.image!(ax, pred; colormap=:navia, colorrange=cr, interpolate=false)
            Colorbar(layout[1, i, Right()], p)
        end

        # Endemics
        extinct_keys = propertynames(frame[].endemic_presence[1])
        n_extinct = length(extinct_keys)
        extinct_strings = collect(string.(extinct_keys))
        menus = map(1:ncols) do i
            Menu(layout[3, i]; default=extinct_strings[i], options=extinct_strings)
        end
        extinct_axes = map(1:ncols, menus) do i, m
            title = lift(m.selection) do s
                _title(s)
            end
            ax = Axis(layout[2, i]; title)
        end
        hidexdecorations!.(extinct_axes; grid=false)
        hideydecorations!.(extinct_axes; grid=false)
        extincts = map(1:ncols) do i
            Observable(replace_missing(rebuild(init.pred_pop, Float32.(getindex.(frame[].endemic_presence, i))), NaN32))
        end

        foreach(extincts, menus) do extinct, menu
            onany(frame, menu.selection) do f, selection
                i = findfirst(==(selection), extinct_strings)
                extinct[] .= replace(Float32.(getindex.(f.endemic_presence, i)), 0.0f0 => NaN32)
                notify(extinct)
            end
        end

        endemic_cmaps = map(1:ncols) do i
            cgrad(ColorScheme([RGB{Float64}(0.0, 0.0, 0.0), RGB{Float64}(i, 0.1i, 1/i)]), 2, categorical=true)
        end
        foreach(extinct_axes, extincts, endemic_cmaps) do ax, extinct, colormap
            Makie.image!(ax, extinct; colorrange=(0.0, 1.0), colormap, interpolate=false)
        end

        # Link
        linkaxes!(pred_axes..., extinct_axes..., ax_lc)
        return nothing
    end
end

function mk_pred(init, ruleset; maxpops=zero(eltype(init.pred_pop)), landcover, tspan, kw...)
    MakieOutput(init;
        kw...,
        tspan,
        fps=100,
        store=false,
        ruleset,
        sim_kw=(; printframe=true),
    ) do (; layout, frame, time)

        colorrange_obs = map(x -> Observable((zero(x), oneunit(x))), maxpops)
        # Landcover
        lc = lift(time) do i
            replace_missing(landcover[Ti(Near(tspan[i]))], NaN32)
        end
        ax_lc = Axis(layout[1, 1]; title="Landcover")
        hidexdecorations!(ax_lc; grid=false)
        hideydecorations!(ax_lc; grid=false)
        Makie.image!(ax_lc, lc; colormap=:batlow, colorrange=(0, 6), interpolate=false)

        # Predators
        pred_keys = propertynames(frame[].pred_pop[1])
        npreds = length(pred_keys)
        ncols = npreds + 1
        pred_axes = map(enumerate(CartesianIndices((2, ncols ÷ 2))[2:end])) do (i, I)
            Axis(layout[Tuple(I)...]; title=_title(pred_keys[i]))
        end
        hidexdecorations!.(pred_axes; grid=false)
        hideydecorations!.(pred_axes; grid=false)
        predators = map(1:npreds) do i
            Observable(rebuild(init.pred_pop, (x -> iszero(x) ? NaN : Float64(x)).(getindex.(frame[].pred_pop, i))))
        end
        on(frame) do f
            foreach(predators, 1:npreds) do  pred, i
                foreach(maximum(f), colorrange_obs) do m, obs 
                    obs[] = m
                    notify(obs)
                end
                pred[] .= (x -> iszero(x) ? NaN : Float64(x)).(getindex.(f.pred_pop, i))
                notify(pred)
            end

        end
        foreach(pred_axes, predators, COLORMAPS[1:npreds], colorrange_obs) do ax, pred, colormap, cr
            Makie.image!(ax, pred; colormap=:navia, colorrange=(zero(mp), cr), interpolate=false)
        end

        # Link
        linkaxes!(ax_lc, pred_axes...)
        return nothing
    end
end

function mk_endemic(init, ruleset; 
    maxpops = zero(eltype(init.pred_pop)),
    landcover, pred_pop, tspan, kw...
)
    MakieOutput(init;
        kw...,
        tspan,
        fps=100,
        store=false,
        ruleset,
        sim_kw=(; printframe=true),
    ) do (; layout, frame, time)
          
        colorrange_obs = map(x -> Observable((zero(x), oneunit(x))), maxpops)
        # Landcover
        lc = lift(time) do i
            replace_missing(landcover[Ti(Near(tspan[i]))], NaN32)
        end
        ax_lc = Axis(layout[1, 1]; title="Landcover")
        # hidexdecorations!(ax_lc; grid=false)
        # hideydecorations!(ax_lc; grid=false)
        Makie.image!(ax_lc, lc; colormap=:batlow, colorrange=(0, 6), interpolate=false)

        # Predators
        pred_keys = propertynames(pred_pop[1])
        npreds = length(pred_keys)
        ncols = npreds + 1
        pred_axes = map(2:ncols) do i
            Axis(layout[1, i]; title=_title(pred_keys[i-1]))
        end
        hidexdecorations!.(pred_axes; grid=false)
        hideydecorations!.(pred_axes; grid=false)
        pred_pop1 = view(pred_pop, Ti=1)
        predators = map(1:npreds) do i
            Observable(rebuild(pred_pop1, (x -> iszero(x) ? NaN : Float64(x)).(getindex.(pred_pop1, i))))
        end
        on(time) do t
            foreach(predators, 1:npreds, colorrange_obs) do pred, i, cr_obs
                pred_pop_t = view(pred_pop, Ti=t)
                pred[] .= (x -> iszero(x) ? NaN : Float64(x)).(getindex.(pred_pop_t, i))
                m = maximum(pred[])
                cr_obs[] = (cr_obs[][1], max(m, cr_obs[][2]))
                notify(cr_obs)
                notify(pred)
            end
        end
        foreach(pred_axes, predators, pred_keys, colorrange_obs) do ax, pred, k, colorrange
            Makie.image!(ax, pred; colormap=:navia, colorrange, interpolate=false)
        end

        # Endemics
        extinct_keys = propertynames(frame[].endemic_presence[1])
        n_extinct = length(extinct_keys)
        extinct_strings = collect(string.(extinct_keys))
        menus = map(1:ncols) do i
            Menu(layout[3, i]; default=extinct_strings[i], options=extinct_strings)
        end
        extinct_axes = map(1:ncols, menus) do i, m
            title = lift(m.selection) do s
                _title(s)
            end
            ax = Axis(layout[2, i]; title)
        end
        hidexdecorations!.(extinct_axes; grid=false)
        hideydecorations!.(extinct_axes; grid=false)
        extincts = map(1:ncols) do i
            Observable(replace_missing(rebuild(init.pred_pop, Float32.(getindex.(frame[].endemic_presence, i))), NaN32))
        end

        foreach(extincts, menus) do extinct, menu
            onany(frame, menu.selection) do f, selection
                i = findfirst(==(selection), extinct_strings)
                extinct[] .= replace(Float32.(getindex.(f.endemic_presence, i)), 0.0f0 => NaN32)
                notify(extinct)
            end
        end

        endemic_cmaps = map(1:ncols) do i
            cgrad(ColorScheme([RGB{Float64}(0.0, 0.0, 0.0), RGB{Float64}(i, 0.1i, 1/i)]), 2, categorical=true)
        end
        foreach(extinct_axes, extincts, endemic_cmaps) do ax, extinct, colormap
            Makie.image!(ax, extinct; colorrange=(0.0, 1.0), colormap, interpolate=false)
        end

        # Link
        linkaxes!(pred_axes..., extinct_axes..., ax_lc)
        return nothing
    end
end

_title(s) = replace(titlecase(string(s)), "_"=>" ")
