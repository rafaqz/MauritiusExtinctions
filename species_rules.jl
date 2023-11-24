using DynamicGrids, Dispersal, LandscapeChange, ModelParameters
using StaticArrays
using Rasters
using ThreadsX
const DG = DynamicGrids
const DD = DimensionalData

# extend(A; to) = _extend(A, to)
# _extend(A, to::Dimension...) = extend(A, (dims...,))
# function _extend(A, to::Tuple)
#     newdims = DD.set(DD.dims(A), map(=>, DD.dims(A, dims), dims)...)
#     newA = similar(A, newdims...)
#     selectors = map(dims) do d
#         rebuild(d, ..(DD.bounds(d)))
#     end
#     newA[selectors...] .= A
#     return newA
# end


# Priors
# Cats prefere prey under 250g (Parsons 2018)

const NV = NamedVector

struct InteractiveGrowth{R,W,GR,CC,I,Su,TS,St} <: Dispersal.GrowthRule{R,W}
    rate::GR
    carrycap::CC
    interactions::I
    suitabilities::Su
    timestep::TS
    nsteps::St
end
function InteractiveGrowth{R,W}(;
    rate,
    carrycap,
    interactions,
    suitabilities,
    timestep,
    nsteps_type=Float64,
) where {R,W}
    InteractiveGrowth{R,W}(
        rate, carrycap, interactions, suitabilities, timestep, zero(nsteps_type)
    )
end

DynamicGrids.modifyrule(rule::InteractiveGrowth, data) = Dispersal.precalc_nsteps(rule, data)

@inline function DynamicGrids.applyrule(data, rule::InteractiveGrowth, populations, I)
    growth_rates = get(data, rule.rate, I...) .* rule.nsteps
    local_suitabilities = map(rule.suitabilities) do (param, val)
        param => get(data, val, I)
    end

    scaling = scale_carrycap(populations, rule.carrycap, rule.interactions, local_suitabilities)
    final_carrycap = rule.carrycap .* scaling

    newpops = map(populations, growth_rates, final_carrycap) do N, rt, k
        if rt > zero(rt)
            if k == zero(k)
                zero(N)
            else
                (N * k) / (N + (k - N) * exp(-rt))
            end
        else
            N * exp(rt)
        end
    end
    if any(isnan, newpops)
        @show populations growth_rates scaling final_carrycap newpops
        error()
    end
    return max.(zero(eltype(populations)), newpops)
end

function scale_carrycap(populations, carrycaps, interactions, suitabilities)
    interaction_factor = sum.(map(i -> i .* populations ./ carrycaps, interactions))
    # @show interaction_factor carrycaps populations

    suitability_factors = map(suitabilities) do (weight, local_intensity)
        weight .* local_intensity
    end
    # @show suitability_factors

    # Simple additive model
    carrycap_adjustment = 1 .+ interaction_factor .+ reduce(.+, suitability_factors)

    # Limit minimum adjustment to 0.0
    return max.(0.0, carrycap_adjustment)
end



function def_syms(
    pred_df, aggfactor, dems, masks, slope_stacks, island_extinct_tables, aux
)
    aggscale = aggfactor^2
    moore = Moore{3}()
    ns_extinct = map(nrow, island_extinct_tables)
    extinct_keys = map(island_extinct_tables) do extinct_table
        Tuple(Symbol.(replace.(extinct_table.Species, Ref(' ' => '_'))))
    end
    ExtinctNVs = map(extinct_keys) do ek
        NamedVector{ek,length(ek)}
    end

    #= Assumptions
    1. cats suppress rodents to some extent, black rats more than norway rats (size selection - norway rats are above 250g)
    2. cats live near people with maybe 2 orders of magnitude higher density than far from people
    3. black rats live anywhere including forest (they're efficient and good climbers), but are outcompeted by norway rats in cities and/or in the presence of cats
    4. norway rats could live everywhere but lose in the forest in competition with black rats, because theyre bad climbers and are less efficient in feeding/metabolism overall. They dominate in coastal areas because theyre better around water, and in high resource conditions because they're bigger.
    5. mice are also implicated in bird deaths, just to make things worse... they outcompete rats in farmland but not forests, and coexist in urban areas.
    6. pigs probably don't care about any of these things and just go wherever they like.
    =#

    # Parameters
    carrycap = NV(
        cat =        0.02,
        black_rat =  30.0,
        norway_rat = 15.0,
        mouse =      52.0,
        pig =        0.3,
        snake =      20.0,
    ) .* aggscale
    interactions = NV(
        cat =        NV(cat= 0.0, black_rat= 0.5, norway_rat= 0.2, mouse= 0.6, pig=0.0, snake= 0.0),
        black_rat =  NV(cat=-0.2, black_rat= 0.0, norway_rat=-0.3, mouse=-0.1, pig=0.0, snake= 0.01),
        norway_rat = NV(cat=-0.1, black_rat=-0.2, norway_rat= 0.0, mouse=-0.1, pig=0.0, snake= 0.01),
        mouse =      NV(cat=-0.2, black_rat=-0.1, norway_rat=-0.1, mouse= 0.0, pig=0.0, snake=-0.1),
        pig =        NV(cat= 0.0, black_rat= 0.0, norway_rat= 0.0, mouse= 0.0, pig=0.0, snake= 0.0),
        snake =      NV(cat=-0.1, black_rat= 0.2, norway_rat= 0.1, mouse= 0.3, pig=0.0, snake= 0.0),
    )
    # 0 is no interaction, less than one negative interaction, more than one positive
    human_dependency = NV(
        cat =        10.0,
        black_rat =   2.0,
        norway_rat =  5.0,
        mouse =       3.0,
        pig =        -0.9,
        snake =      -0.9,
    )
    forest_preference = Param.(NV(
        cat =        -0.2,
        black_rat =   0.4,
        norway_rat = -0.2,
        mouse =      -0.2,
        pig =         0.5,
        snake =       0.3,
    ))
    spread_rate = NV(
        cat =         10.0,
        black_rat =   1.0,
        norway_rat =  1.0,
        mouse =       0.5,
        pig =         10.0,
        snake =       1.0,
    )

    # Example data
    human_intensity = 0.0
    forest_density = 0.9

    # How much these species are supported outside of this system
    # populations = NV(cat=0.01, black_rat=25.0, norway_rat=10.0, pig=0.3)

    # suitabilities = (human_dependency => human_intensity, forest_preference => forest_density)
    # scale_carrycap(populations, carrycap, interactions, suitabilities)

    # TODO :use forest density here
    suitabilities = (forest_preference => Aux{:never_cleared}(),)

    ag_masks = map(masks) do mask
        Rasters.aggregate(Rasters.Center(), mask, aggfactor)
    end
    ag_slope = map(slope_stacks) do s
        Rasters.aggregate(maximum, replace_missing(s.slope, 0), aggfactor)
    end
    A = first(aux)
    map(aux) do A
        ag_A1 = DimensionalData.modify(BitArray, Rasters.aggregate(Center(), A, (X(aggfactor), Y(aggfactor))))
        # Extend a century earlier
        ag_A = Rasters.extend(maybeshiftlocus(ag_A1; to=(Ti(1500:1:2020),))
        broadcast_dims!(identity, view(ag_A, Ti=1500..1600), view(ag_A1, Ti=At(1600)))
    end

    ag_roughness = map(dems) do dem
        # Clip roughness at a maximum of 500
        rgh = min.(roughness(replace_missing(dem, 0)), 500)
        Rasters.aggregate(maximum, rgh ./ 500, aggfactor)
    end
    ag_dems = map(dems) do dem
        Rasters.aggregate(Rasters.Center(), dem, aggfactor)
    end
    ag_local_cost = map(ag_dems) do dem
        window = Window{2}()
    end

    kernel = DispersalKernel(
        stencil=moore,
        formulation=ExponentialKernel()
    )
    stencil_masks = map(ag_masks) do mask
        StencilArray(mask, kernel; padding=Halo{:out}())
    end
    spreadability = map(ag_roughness, ag_slope, ag_masks) do r, s, m
        # A = (1 .- s) .* m
        A = (1 .- sqrt.(s)) .*  m
        StencilArray(A, moore; padding=Halo{:out}())
    end
    stencil_dems = map(dems) do dem
        A = Rasters.aggregate(Rasters.Center(), replace_missing(dem, 0), aggfactor)
        StencilArray(A, moore; padding=Halo{:out}())
    end

    pred_keys = Tuple(Symbol.(pred_df.name))
    PredNV = NamedVector{pred_keys,length(pred_keys)}
    pred_names = PredNV(pred_keys)
    pred_indices = PredNV(ntuple(identity, length(pred_keys)))
    island_names = NamedTuple{keys(masks)}(keys(masks))

    one_individual = map(pred_indices) do _
        1.0
    end

    introductions = map(island_names) do key
        years = pred_df[!, "$(key)_year"]
        lons = pred_df[!, "$(key)_lon"]
        lats = pred_df[!, "$(key)_lat"]
        inits = map(pred_indices, one_individual) do i, oi
            x = zeros(length(pred_keys))
            x[i] = oi * 50 # Lets just say ten of everything where introduced...
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

    pred_growth = InteractiveGrowth{:pred_pops}(;
        rate=pred_rmax,
        carrycap,
        interactions,
        suitabilities,
        timestep=1,
    )

    habitat_requirements = map(ExtinctNVs) do ExtinctNV
        rand(ExtinctNV)
    end
    hunting_suscepts = map(ExtinctNVs) do ExtinctNV
        rand(ExtinctNV)
    end
    # Dummy values. This will be calculated by the optimiser
    pred_suscepts = map(ns_extinct, ExtinctNVs) do n_extinct, ExtinctNV
        cat_sus = ExtinctNV(LinRange(0.0, 0.02, n_extinct))
        black_rat_sus = ExtinctNV(LinRange(0.0, 0.01, n_extinct))
        norway_rat_sus = ExtinctNV(LinRange(0.0, 0.01, n_extinct))
        mouse_sus = ExtinctNV(LinRange(0.0, 0.001, n_extinct))
        pig_sus = ExtinctNV(LinRange(0.0, 0.005, n_extinct))
        snake_sus = ExtinctNV(LinRange(0.0, 0.005, n_extinct))
        pred_suscept = map(cat_sus, black_rat_sus, norway_rat_sus, mouse_sus, pig_sus, snake_sus) do c, br, nr, m, p, s
            PredNV((c, br, nr, m, p, s)) ./ aggscale
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
                @show intro.year current_year
                p = intro.geometry
                I = DimensionalData.dims2indices(D, (X(Contains(p.X)), Y(Contains(p.Y))))
                pred_pops[I...] = pred_pops[I...] .+ intro.init
                @show pred_pops[I...]
            end
        end
    end

    risks = let stochastic_extirpation=0.001*aggfactor, f=tanh
        Cell{Tuple{:pred_pops,:endemic_presences},:endemic_presences}() do data, (pred_pops, endemic_presences), I
            pred_suscept = DG.aux(data).pred_suscept
            map(endemic_presences, pred_suscept) do present, ps
                if present
                    # rand() < hp * hs & rand() < sum(map(*, ps, pred_pops))
                    rand() > f(sum(map(*, ps, pred_pops))) && rand() > stochastic_extirpation
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
            hp = get(data, Aux{:never_cleared}(), I)
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
            uc = get(data, Aux{:never_cleared}(), I)
            presences .& uc
        end
    end

    pred_kernels = map(spread_rate) do s
        DispersalKernel(
            stencil=moore,
            # formulation=ExponentialKernel(Param(s, bounds=(0.000000000001, 100.0)))
            formulation=ExponentialKernel(s),
        )
    end
    pred_spread = let aggfactor=aggfactor, pred_kernels=pred_kernels
        SetNeighbors{:pred_pops}(pred_kernels[1]) do data, hood, Ns, I
            any(isnan, Ns) && error("pred_spread N is: $Ns")
            Ns == zero(Ns) && return nothing
            dem = DG.aux(data)[:dem]
            spr = DG.mask(data)#[:spreadability]
            spr_nbrs = DG.neighbors(spr, I)
            elev_nbrs = DG.neighbors(dem, I)
            @inbounds elev_center = dem[I...]

            # spreadability = DG.neighbors(DG.mask(data), I)
            sum = zero(Ns) # TODO needs cell size here
            cellsize = (30 * aggfactor)
            any(isnan, spr_nbrs) && return zero(Ns)

            start = rand(0:length(hood)-1)
            @inbounds for ix in eachindex(hood)
                i = start + ix
                if i > length(hood)
                    i = i - length(hood)
                end
                Ih = DG.indices(hood, I)[i]
                ks = getindex.(DynamicGrids.kernel.(pred_kernels), i)
                d = DG.distances(hood)[i]
                sp = spr_nbrs[i]
                e = elev_nbrs[i]
                # any(isnan, N) && error("N is NaN")
                # isnan(k) && error("k is NaN")
                # isnan(sp) && error("sp is NaN")
                # slope = (e - elev_center) / (cellsize * d)
                # Imhof Tobler equation
                # speed = (6ℯ.^(slopefactor .* (slope .+ distfactor)))

                # slope_factor = 1 - min(slope, 1) # TODO needs cell size here
                slope_factor = 1.0000000001 - min(abs(e - elev_center) / cellsize, 1) # TODO needs cell size here
                # @show N  k  sp slope_factor one_individual
                # propagules = N .* k .* sp .* slope_factor .+ ((rand(typeof(N)) .- 0.5) .* one_individual)
                # propagules = Ns * sp * slope_factor .* kr .* rand(typeof(Ns))
                propagules = trunc.(Ns * sp * slope_factor .* rand(typeof(Ns)) .^ 3 .* ks .* 40)
                any(isnan, propagules) && error("NaNs found Ns: $Ns sp: $sp slope_factor: $slope_factor ks: $ks")
                sum1 = sum + propagules
                # If we run out of propagules
                if any(sum1 .> Ns)
                    propagules = min.(sum1, Ns) .- sum
                    sum = sum + propagules
                else
                    sum = sum1
                end
                add!(data[:pred_pops], propagules, Ih...)
            end
            @inbounds any(isnan, data[:pred_pops][I...]) && error("nan in output")
            @inbounds sub!(data[:pred_pops], sum, I...)
            return nothing
        end
    end # let
    pred_spread

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
        1549:2018
    end
    inits = map(pred_pops, endemic_presences) do pred_pops, endemic_presences
        (; pred_pops, endemic_presences)
    end
    ruleset = Ruleset(
        introduction_rule,
        pred_spread, pred_growth,
        endemic_recouperation, risks, clearing;
        boundary=Use()
    )

    pred_ruleset = Ruleset(
        introduction_rule,
        pred_spread,
        pred_growth;
        boundary=Use()
    )

    outputs_kw = map(island_names, ns_extinct, tspans) do island, n_extinct, tspan
        aux = (;
            spreadability=getproperty(spreadability, island),
            introductions=getproperty(introductions, island),
            pred_suscept=getproperty(pred_suscepts, island),
            recouperation_rate=getproperty(recouperation_rates, island),
            habitat_requirement=getproperty(habitat_requirements, island),
            dem=getproperty(stencil_dems, island),
            never_cleared=ag_never_cleared,
            forested=ag_forested,
            urban=ag_urban,
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
    pred_inits = map(pred_pops) do pred_pops
        (; pred_pops)
    end
    pred_outputs = map(pred_inits, outputs_kw) do init, kw
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

    return (; ruleset, pred_ruleset, inits, pred_inits, outputs, pred_outputs, outputs_kw, ag_masks)#, last_obs
end

function mk(init, ruleset; kw...)
    MakieOutput(init;
        kw...,
        fps=100,
        store=false,
        ruleset,
        sim_kw=(; printframe=true),
    ) do fig, frame, time
        pred_keys = propertynames(frame[].pred_pops[1])
        ncols = length(pred_keys)
        extinct_keys = propertynames(frame[].endemic_presences[1])
        n_extinct = length(extinct_keys)
        sg = SliderGrid(fig[3,1:2],
            map(1:ncols) do i
                (startvalue=i, range=1:n_extinct, format=i -> string(extinct_keys[i]))
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
        foreach(pred_axes, predators, pred_keys) do ax, pred, k
            Makie.image!(ax, pred; colormap=:magma, interpolate=false)
        end
        foreach(extinct_axes, extincts, Iterators.Cycle([:blues, :reds, :greens]), sg.labels) do ax, extinct, colormap, label
            Makie.image!(ax, extinct; colorrange=(0.0, 1.0), colormap, interpolate=false)
        end
        return nothing
    end
end

function mk_pred(init, ruleset, mask; kw...)
    MakieOutput(init;
        kw...,
        fps=100,
        store=false,
        ruleset,
        sim_kw=(; printframe=true),
    ) do fig, frame, time
        # nanmask = map(x -> x ? 1.0f0 : NaN32, mask)
        pred_keys = propertynames(frame[].pred_pops[1])
        ncols = length(pred_keys)
        pred_axes = map(1:ncols) do i
            Axis(fig[1, i]; title=string(pred_keys[i]))
        end
        linkaxes!(pred_axes...)
        predators = map(1:ncols) do i
            Observable(Array(Float32.(getindex.(frame[].pred_pops, i))))
        end

        on(frame) do f
            foreach(predators, 1:ncols) do pred, i
                pred[] .= getindex.(f.pred_pops, i)
                notify(pred)
            end
        end
        foreach(pred_axes, predators, pred_keys) do ax, pred, k
            Makie.image!(ax, pred; colormap=:viridis, interpolate=false)
        end
        return nothing
    end
end


