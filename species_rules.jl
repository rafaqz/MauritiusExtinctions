using DynamicGrids, Dispersal, LandscapeChange, ModelParameters
using StaticArrays
using Rasters
using ThreadsX
const DG = DynamicGrids

# Priors
# Cats prefere prey under 250g (Parsons 2018)

function def_syms(
    pred_df, aggfactor, dems, masks, slope_stacks, island_extinct_tables
)
    window = Window{2}()
    ns_extinct = map(nrow, island_extinct_tables)
    extinct_keys = map(island_extinct_tables) do extinct_table
        Tuple(Symbol.(replace.(extinct_table.Species, Ref(' ' => '_'))))
    end
    ExtinctNVs = map(extinct_keys) do ek
        NamedVector{ek,length(ek)}
    end

    ag_masks = map(masks) do mask
        Rasters.aggregate(Rasters.Center(), mask, aggfactor)
    end
    ag_slope = map(slope_stacks) do s
        Rasters.aggregate(maximum, replace_missing(s.slope, 0), aggfactor)
    end
    ag_uncleared = modify(BitArray, Rasters.aggregate(Center(), uncleared, (X(aggfactor), Y(aggfactor))))

    ag_roughness = map(dems) do dem
        # Clip roughness at a maximum of 500
        rgh = min.(roughness(replace_missing(dem, 0)), 500)
        Rasters.aggregate(maximum, rgh ./ 500, aggfactor)
    end

    one_individual = NamedVector(
        cats=0.0001 * aggfactor, 
        black_rats=0.00001 * aggfactor, 
        pigs=0.0001 * aggfactor,
    )
    kernel = DispersalKernel(
        stencil=window,
        formulation=ExponentialKernel()
    )
    stencil_masks = map(ag_masks) do mask
        StencilArray(mask, kernel; padding=Halo{:out}())
    end
    spreadability = map(ag_roughness, ag_slope, ag_masks) do r, s, m
        # A = (1 .- s) .*  m
        A = (1 .- sqrt.(s)) .*  m
        StencilArray(A, window; padding=Halo{:out}())
    end
    stencil_dems = map(dems) do dem
        A = Rasters.aggregate(Rasters.Center(), replace_missing(dem, 0), aggfactor)
        StencilArray(A, window; padding=Halo{:out}())
    end

    pred_keys = Tuple(Symbol.(pred_df.name))
    PredNV = NamedVector{pred_keys,length(pred_keys)}
    pred_names = PredNV(pred_keys)
    pred_indices = PredNV(ntuple(identity, length(pred_keys)))
    island_names = NamedTuple{keys(masks)}(keys(masks))

    introductions = map(island_names) do key
        years = pred_df[!, "$(key)_year"]
        lons = pred_df[!, "$(key)_lon"]
        lats = pred_df[!, "$(key)_lat"]
        inits = map(pred_indices, one_individual) do i, oi
            x = zeros(length(pred_keys))
            x[i] = oi * 10 # Lets just say ten of everything where introduced...
            PredNV(x)
        end
        map(years, lons, lats, inits) do year, lat, lon, init
            (; year, geometry=(X=lon, Y=lat), init)
        end |> Tuple |> NamedTuple{pred_keys}
    end

    pred_rmax = PredNV(pred_df.rmax)
    pred_max_density = PredNV(pred_df.max_density)
    pred_pops = map(ag_masks) do m
        map(_ -> map(_ -> 0.0, pred_rmax), m)
    end

    habitat_requirements = map(ExtinctNVs) do ExtinctNV
        rand(ExtinctNV)
    end
    hunting_suscepts = map(ExtinctNVs) do ExtinctNV
        rand(ExtinctNV)
    end
    # Dummy values. This will be calculated by the optimiser
    pred_suscepts = map(ns_extinct, ExtinctNVs) do n_extinct, ExtinctNV
        catsus = ExtinctNV(LinRange(0.0, 0.03, n_extinct))
        ratsus = ExtinctNV(LinRange(0.0, 0.02, n_extinct))
        pigsus = ExtinctNV(LinRange(0.0, 0.01, n_extinct))
        pred_suscept = map(catsus, ratsus, pigsus) do c, r, p
            PredNV((c, r, p))
        end
    end

    # Every species is everywhere initially, in this dumb model
    endemic_presences = map(ag_masks, hunting_suscepts) do mask, hunting_suscept
        map(mask) do m
            map(_ -> m, hunting_suscept) 
        end
    end

    #### Rules ##########################################################333

    introduction_rule = SetGrid{:pred_pops}() do data, pred_pops, t
        D = dims(DG.init(data).pred_pops)
        current_year = currenttime(data)
        intros = DG.aux(data)[:introductions]
        foreach(intros) do intro
            if intro.year == current_year
                p = intro.geometry
                I = DimensionalData.dims2indices(D, (X(Contains(p.X)), Y(Contains(p.Y))))
                pred_pops[I...] = pred_pops[I...] .+ intro.init 
            end
        end
    end

    risks = let stochastic_extirpation=0.001*aggfactor
        Cell{Tuple{:pred_pops,:endemic_presences},:endemic_presences}() do data, (pred_pops, endemic_presences), I
            pred_suscept = DG.aux(data).pred_suscept
            map(endemic_presences, pred_suscept) do present, ps
                if present
                    # rand() < hp * hs & rand() < sum(map(*, ps, pred_pops))
                    rand() > sum(map(*, ps, pred_pops)) &&
                    rand() > stochastic_extirpation
                    # ... etc 
                else
                    false
                end
            end
        end
    end
    # hp = get(data, Aux{:hunting_pressure}(), I)

    habitat = let
        Cell{:presences}() do data, presences, I
            hp = get(data, Aux{:uncleared}(), I)
            habitat_requirement = DG.aux(data).habitat_requirement
            map(presences, habitat_requirement) do present, hs
                if present
                    rand() < hp * hs
                else
                    false
                end
            end
        end
    end

    clearing = let
        Cell{:endemic_presences}() do data, presences, I
            uncleared = get(data, Aux{:uncleared}(), I)
            map(presences) do present
                present & uncleared
            end
        end
    end

    pred_spread = let aggfactor=aggfactor, one_individual=one_individual
        SetNeighbors{:pred_pops}(kernel) do data, hood, N, I 
            N == zero(N) && return nothing
            spr = DG.aux(data)[:spreadability]
            dem = DG.aux(data)[:dem]
            spr_nbrs = DG.neighbors(spr, I)
            dem_nbrs = DG.neighbors(dem, I)
            dem_center = dem[I...]

            # spreadability = DG.neighbors(DG.mask(data), I)
            sum = zero(N)
            for (i, k, sp, d) in zip(DG.indices(hood, I), DG.kernel(hood), spr_nbrs, dem_nbrs)
                # any(isnan, N) && error("N is NaN")
                # isnan(k) && error("k is NaN")
                # isnan(sp) && error("sp is NaN")
                slope_factor = 1 - min(abs(d - dem_center) / (15 * aggfactor), 1) # TODO needs cell size here
                # @show N  k  sp slope_factor one_individual
                propagules = N .* k .* sp .* slope_factor .+ ((rand(typeof(N)) .- 0.5) .* 100.0 .* one_individual)
                # any(isnan, propagules) && error()
                @inbounds add!(data[:pred_pops], propagules, i...)
                sum += propagules
            end
            @inbounds sub!(data[:pred_pops], sum, I...)
            return nothing
        end
    end # let

    pred_growth = LogisticGrowth{:pred_pops}(; 
        rate=pred_rmax,
        carrycap=pred_max_density,
        timestep=1.0,
    )
    # This is not really allee extinction, just assuring
    # at least one individual exists.
    pred_allee = AlleeExtinction{:pred_pops}(;
        minfounders=one_individual,
    )

    recouperation_rates = map(ExtinctNVs) do ExtinctNV
        rand(ExtinctNV) / 10
    end
    endemic_recouperation = let 
        Neighbors{:endemic_presences}(Moore(2)) do data, hood, prescences, I
            any(prescences) || return prescences
            recouperation_rate = DG.aux(data).recouperation_rate
            i = 0
            nbr_sums = sum(hood)
            map(prescences, nbr_sums, recouperation_rate) do p, n_nbrs, rr
                if p 
                    true
                else
                    rand() < (n_nbrs / length(hood) * rr)
                end
            end
        end
    end

    tspans = map(introductions) do intros
        # t1 = minimum(x -> x.year, intros) - 2
        1600:2018
    end
    inits = map(pred_pops, endemic_presences) do pred_pops, endemic_presences
        (; pred_pops, endemic_presences)
    end
    ruleset = Ruleset(
        introduction_rule, 
        pred_spread, 
        pred_allee, pred_growth, 
        endemic_recouperation, risks, clearing; 
        boundary=Ignore()
    )
    # ruleset = Ruleset(Chain(endemic_recouperation, risks, clearing); boundary=Ignore())

    outputs_kw = map(island_names, ns_extinct, tspans) do island, n_extinct, tspan
        aux = (; 
            spreadability=getproperty(spreadability, island), 
            introductions=getproperty(introductions, island),
            pred_suscept=getproperty(pred_suscepts, island),
            recouperation_rate=getproperty(recouperation_rates, island),
            habitat_requirement=getproperty(habitat_requirements, island),
            dem=getproperty(stencil_dems, island),
            uncleared=ag_uncleared,
        )
        (; aux, mask=getproperty(stencil_masks, island), tspan, n_extinct)
    end
    outputs = map(inits, outputs_kw) do init, kw
        ArrayOutput(init; kw...) 
        # trans_output = TransformedOutput(init; kw...) do f
        #     Base.reduce(parent(parent(f.endemic_presences)); init=zero(eltype(f.endemic_presences))) do acc, xs 
        #         map(|, acc, xs)
        #     end
        # end
    end
    # Sum as we go
    # ncells = sum(DG.mask(output))
    # last_obs = map(island_names) do island
    #     ExtinctNV(map(first, getproperty(island_extinct_tables, island).extinct))
    # end

    return ruleset, inits, outputs, outputs_kw#, last_obs
end

function mk(init, ruleset; kw...)
    MakieOutput(init;
        kw...,
        fps=100, 
        store=false,
        ruleset,
        sim_kw=(; printframe=true),
    ) do fig, time, frame
        pred_keys = propertynames(frame[].pred_pops[1])
        ncols = length(pred_keys)
        extinct_keys = propertynames(frame[].endemic_presences[1])
        n_extinct = length(extinct_keys) 
        sg = SliderGrid(fig[3,1:2], 
            map(1:ncols) do i
                (startvalue=i, range=1:n_extinct, format=i ->string(extinct_keys[i]))
            end...
        )
        pred_axes = map(1:ncols) do i
            Axis(fig[1, i]; title=string(pred_keys[i]))
        end
        extinct_axes = map(1:ncols, sg.labels) do i, l
            ax = Axis(fig[2, i]; title=l.text)
        end
        linkaxes!(pred_axes..., extinct_axes...)
        predators = map(1:ncols) do i
            Observable(Array(getindex.(frame[].pred_pops, i)))
        end
        extincts = map(1:ncols) do i
            Observable(Array(getindex.(frame[].endemic_presences, i)))
        end
       
        on(frame) do f
            foreach(predators, 1:ncols) do  pred, i
                pred[] .= getindex.(f.pred_pops, i)
                notify(pred)
            end
        end
        foreach(extincts, sg.sliders) do extinct, slider
            onany(frame, slider.value) do f, s
                extinct[] .= getindex.(f.endemic_presences, s)
                notify(extinct)
            end
        end
        foreach(pred_axes, predators, [:magma, :viridis, :cividis], pred_keys) do ax, pred, colormap, k
            Makie.image!(ax, pred; colorrange=(0.0, 1.0), colormap, interpolate=false)
        end
        foreach(extinct_axes, extincts, [:blues, :reds, :greens], sg.labels) do ax, extinct, colormap, label
            Makie.image!(ax, extinct; colorrange=(0.0, 1.0), colormap, interpolate=false)
        end
        return nothing
    end
end


