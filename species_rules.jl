using DynamicGrids, Dispersal, LandscapeChange, ModelParameters
using StaticArrays
using Rasters
using ThreadsX
const DG = DynamicGrids
const DD = DimensionalData

# Priors
# Cats prefere prey under 250g (Parsons 2018)

const NV = NamedVector

struct InteractiveCarryCap{R,W,CC,CS,I} <: Dispersal.GrowthRule{R,W}
    carrycap::CC
    carrycap_scaling::CS
    inputs::I
end
function InteractiveCarryCap{R,W}(; carrycap, carrycap_scaling, inputs=(;),) where {R,W}
    InteractiveCarryCap{R,W}(carrycap, carrycap_scaling, inputs)
end

Base.@assume_effects :foldable function DynamicGrids.applyrule(
    data, rule::InteractiveCarryCap,
    populations::NamedVector{Keys}, I,
) where Keys
    local_inputs = get(data, rule.inputs, I)
    relative_pop = NamedTuple(populations ./ rule.carrycap)
    # combine populations
    params = merge(relative_pop, local_inputs)

    scaling = map(rule.carrycap_scaling) do val_f
        f = DynamicGrids._unwrap(val_f)
        (oneunit(eltype(populations)) + f(params)::Float32)::Float32
    end |> NamedVector
    # # Carrycap cant be zero
    new_carrycaps = max.(oneunit(eltype(rule.carrycap)) .* 1f-10, rule.carrycap .* scaling)
    return new_carrycaps
end

struct ExtirpationRisks{R,W,F,P,PS,SE} <: DynamicGrids.CellRule{R,W}
    f::F
    pred_pop::P
    pred_suscept::PS
    stochastic_extirpation::SE
end
function ExtirpationRisks{R,W}(; f, pred_pop, pred_suscept, stochastic_extirpation) where {R,W}
    ExtirpationRisks{R,W}(f, pred_pop, pred_suscept, stochastic_extirpation)
end

@inline function DynamicGrids.applyrule(data, rule::ExtirpationRisks, endemic_presence, I)
    pred_pop = get(data, rule.pred_pop, I)
    pred_suscept = get(data, rule.pred_suscept)
    map(endemic_presence, pred_suscept) do present, ps
        if present
            # rand() < hp * hs & rand() < sum(map(*, ps, pred_pop))
            rand(Float32) > rule.f(sum(map(*, ps, pred_pop))) && rand(Float32) > rule.stochastic_extirpation
            # ... etc
        else
            false
        end
    end
end

function agg_aux(masks::NamedTuple, slope_stacks, dems, lcs, aggfactor)
    map(masks, slope_stacks, dems, lcs) do args...
        agg_aux(args..., aggfactor)
    end
end
function agg_aux(mask_orig::AbstractArray, slope_orig, dem_orig, lc_orig, aggfactor)
    mask = Rasters.aggregate(Rasters.Center(), mask_orig, aggfactor)
    slope = Rasters.aggregate(maximum, replace_missing(slope_orig.slope, 0), aggfactor)
    # Clip roughness at a maximum of 500
    # rgh = Float32.(replace_missing(Rasters.aggregate(maximum, min.(roughness(replace_missing(dem_orig, 0)), 500)./ 500, aggfactor), NaN))
    dem = Float32.(replace_missing(Rasters.aggregate(Rasters.Center(), dem_orig, aggfactor), 0))
    lc_ag = Rasters.aggregate(mean, rebuild(lc_orig; missingval=nothing), (X(aggfactor), Y(aggfactor)))
    lc_ag1 = Rasters.extend(lc_ag; to=(Ti(Sampled(1500:1:2020; sampling=Intervals(Start())))), missingval=false)
    map(DimensionalData.layers(lc_ag1)) do A
        broadcast_dims!(identity, view(A, Ti=1500..1600), view(A, Ti=At(1600)))
    end
    lc = map(CartesianIndices(first(lc_ag1))) do I
        Float32.(NV(lc_ag1[I]))
    end
    (; mask, slope, dem, lc)
end

function def_syms(
    pred_df, introductions_df, island_extinct_tables, auxs, aggfactor;
    replicates=nothing, pred_pops_aux=map(_ -> nothing, dems),
)
    aggscale = aggfactor^2
    window = Window{2}()
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
    carrycap = Float32.(NV(
        cat =        0.02,
        black_rat =  30.0,
        norway_rat = 15.0,
        mouse =      52.0,
        pig =        0.3,
        snake =      20.0,
        macaque =    1.0,
    )) .* aggscale

    pred_funcs = (;
        cat =        p -> 1.0f0p.black_rat + 0.3f0p.norway_rat + 1.0f0p.mouse + 10f0p.urban + 2f0p.cleared,
        black_rat  = p -> -0.2f0p.cat - 0.1f0p.norway_rat - 0.1f0p.mouse + 0.5f0p.native + 0.3f0p.abandoned + 1p.urban,
        norway_rat = p -> -0.1f0p.cat - 0.1f0p.black_rat - 0.1f0p.mouse + 1.5f0p.urban - 0.2f0p.native,
        mouse =      p -> -0.3f0p.cat - 0.2f0p.black_rat - 0.2f0p.norway_rat + 0.8f0p.cleared + 1.5f0p.urban,
        pig =        p -> 0.5f0p.native + 0.4f0p.abandoned - 2f0p.urban - 1.0f0p.cleared,
        snake =      p -> -0.2f0p.cat + 0.2f0p.black_rat + 0.3f0p.mouse - 0.5f0p.urban + 0.3f0p.native,
        macaque =    p -> 1.0f0p.abandoned + 0.7f0p.forestry + 0.4f0p.native - 1.0f0p.urban - 0.8f0p.cleared
    )

    # These need to somewhat balance low growth rates
    spread_rate = NV(
        cat =         30.0f0,
        black_rat =   1.0f0,
        norway_rat =  1.0f0,
        mouse =       0.5f0,
        pig =         15.0f0,
        snake =       1.0f0,
        macaque =     5.0f0,
    )

    # How much these species are supported outside of this system
    # populations = NV(cat=0.01, black_rat=25.0, norway_rat=10.0, pig=0.3)

    # suitabilities = (human_dependency => human_intensity, forest_preference => forest_density)
    # scale_carrycap(populations, carrycap, interactions, suitabilities)

    kernel = DispersalKernel(
        stencil=moore,
        formulation=ExponentialKernel(Param(1.0f0, bounds=(0.0f0, 2.0f0))),
        cellsize=1.0f0,
    )
    stencil_masks = map(auxs) do aux
        StencilArray(aux.mask, kernel; padding=Halo{:out}())
    end
    stencil_dems = map(auxs) do aux
        StencilArray(aux.dem, moore; padding=Halo{:out}())
    end

    pred_keys = Tuple(Symbol.(pred_df.name))
    PredNV = NamedVector{pred_keys,length(pred_keys)}
    pred_names = PredNV(pred_keys)
    pred_indices = PredNV(ntuple(identity, length(pred_keys)))
    island_names = NamedTuple{keys(auxs)}(keys(auxs))
    pred_inits = map(pred_indices) do i
        x = zeros(Float16, length(pred_keys))
        x[i] = 50 * aggfactor
        PredNV(x)
    end

    introductions = map(island_names) do key
        island_df = filter(r -> r.island == string(key), introductions_df)
        map(eachrow(island_df)) do r
            init = pred_inits[Symbol(r.species)]
            (; year=r.year, geometry=(X=r.lon, Y=r.lat), init)
        end
    end

    pred_rmax = Float32.(PredNV(pred_df.rmax))
    pred_max_density = Float32.(PredNV(pred_df.max_density))
    pred_pops = map(auxs) do aux
        map(_ -> map(_ -> 0.0f0, pred_rmax), aux.mask)
    end
    pred_carrycaps = map(auxs) do aux
        map(_ -> carrycap, aux.mask)
    end


    pred_carrycap_rule = InteractiveCarryCap{:pred_pop,:pred_carrycap}(;
        carrycap,
        carrycap_scaling=map(Val, pred_funcs),
        inputs=Aux{:landcover}(),
    )

    pred_growth_rule = LogisticGrowth{:pred_pop}(;
        rate=pred_rmax,
        carrycap=Grid{:pred_carrycap}(),
        timestep=1,
        nsteps_type=Float32,
    )

    habitat_requirements = map(ExtinctNVs) do ExtinctNV
        Float32.(rand(ExtinctNV))
    end
    hunting_suscepts = map(ExtinctNVs) do ExtinctNV
        Float32.(rand(ExtinctNV))
    end
    # Dummy values. This will be calculated by the optimiser
    pred_suscepts = map(ns_extinct, ExtinctNVs) do n_extinct, ExtinctNV
        cat_sus = ExtinctNV(LinRange(Float32(0.0), Float32(0.01), n_extinct))
        black_rat_sus = ExtinctNV(LinRange(Float32(0.0), Float32(0.005), n_extinct))
        norway_rat_sus = ExtinctNV(LinRange(Float32(0.0), Float32(0.07), n_extinct))
        mouse_sus = ExtinctNV(LinRange(Float32(0.0), Float32(0.001), n_extinct))
        pig_sus = ExtinctNV(LinRange(Float32(0.0), Float32(0.005), n_extinct))
        snake_sus = ExtinctNV(LinRange(Float32(0.0), Float32(0.005), n_extinct))
        macaque_sus = ExtinctNV(LinRange(Float32(0.0), Float32(0.003), n_extinct))
        pred_suscept = map(cat_sus, black_rat_sus, norway_rat_sus, mouse_sus, pig_sus, snake_sus, macaque_sus) do c, br, nr, m, p, s, ma
            PredNV((c, br, nr, m, p, s, ma)) ./ aggscale
        end |> collect
    end

    # Every species is everywhere initially, in this dumb model
    endemic_presences = map(auxs, hunting_suscepts) do aux, hunting_suscept
        map(aux.mask) do m
            map(_ -> m, hunting_suscept)
        end
    end


    #### Rules ##########################################################333

    introduction_rule = let introductions_aux=Aux{:introductions}()
        SetGrid{:pred_pop}() do data, pred_pop, t
            D = dims(DG.init(data).pred_pop)
            current_year = currenttime(data)
            intros = get(data, introductions_aux)
            foreach(intros) do intro
                if intro.year == current_year
                    p = intro.geometry
                    I = DimensionalData.dims2indices(D, (X(Contains(p.X)), Y(Contains(p.Y))))[1:2]
                    # x = view(pred_pop, I..., :) .+ (intro.init,)
                    # @show x I size(intro.init) size(pred_pop)
                    pred_pop[I..., :] .= view(pred_pop, I..., :) .+ (intro.init,)
                end
            end
        end
    end

    habitat = let native=Aux{:native}()
        Cell{:presences}() do data, presences, I
            hp = get(data, native, I)
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

    clearing_rule = let landcover=Aux{:landcover}()
        Cell{:endemic_presence}() do data, presences, I
            lc = get(data, landcover, I)
            presences .& (lc.native > 0.2f0)
        end
    end

    pred_kernels = map(spread_rate) do s
        DispersalKernel(
            stencil=moore,
            # formulation=ExponentialKernel(Param(s, bounds=(0.000000000001, 100.0)))
            formulation=ExponentialKernel(s),
            cellsize=1.0f0,
        )
    end
    pred_spread_rule = let demaux=Aux{:dem}(), aggfactor=aggfactor, pred_kernels=pred_kernels, carrycap=carrycap, slope_scalar=Param(5.0, bounds=(1.0, 20.0))
        SetNeighbors{Tuple{:pred_pop,:pred_carrycap}}(pred_kernels[1]) do data, hood, (Ns, _), I
            Ns === zero(Ns) && return nothing
            dem = get(data, demaux)
            reps = DynamicGrids.replicates(data)
            carrycap_nbrs = if isnothing(reps)
                DG.neighbors(DG.grids(data).pred_carrycap, I)
            else
                DG.neighbors(DG.grids(data).pred_carrycap, (I..., reps))
            end
            elev_nbrs = DG.neighbors(dem, I)
            @inbounds elev_center = dem[I...]

            sum = zero(Ns) # TODO needs cell size here
            cellsize = (100 * aggfactor)

            # Randomise hood starting position to avoid directional artifacts in output
            start = rand(0:length(hood)-1)
            @inbounds for ix in eachindex(hood)
                # Rotate indices in relation to starting point
                i = start + ix
                if i > length(hood)
                    i = i - length(hood)
                end
                sp = carrycap_nbrs[i] ./ carrycap
                any(map(isnan, sp)) && error("sp is NaN")
                Ih = DG.indices(hood, I)[i]
                ks = getindex.(DynamicGrids.kernel.(pred_kernels), i)
                d = DG.distances(hood)[i]
                e = elev_nbrs[i]
                # slope = (e - elev_center) / (cellsize * d)
                # Imhof Tobler equation
                # speed = (6ℯ.^(slopefactor .* (slope .+ distfactor)))

                # slope_factor = 1 - min(slope, 1) # TODO needs cell size here
                slope_factor = max(0, 1 - slope_scalar * (elev_center - e) / cellsize)
                # @show N  k  sp slope_factor one_individual
                # propagules = N .* k .* sp .* slope_factor .+ ((rand(typeof(N)) .- 0.5) .* one_individual)
                # propagules = Ns * sp * slope_factor .* kr .* rand(typeof(Ns))
                propagules = trunc.(Float32.(Ns .* sp .* slope_factor .* rand(typeof(Ns)) .^ 3 .* ks .* 40))
                sum1 = sum + propagules
                # If we run out of propagules
                if any(sum1 .> Ns)
                    propagules = min.(sum1, Ns) .- sum
                    sum = sum + propagules
                else
                    sum = sum1
                end
                add!(data[:pred_pop], propagules, Ih...)
            end
            @inbounds sub!(data[:pred_pop], sum, I...)
            return nothing
        end
    end # let

    recouperation_rates = map(ExtinctNVs) do ExtinctNV
        Float32.(rand(ExtinctNV) / 10)
    end
    endemic_recouperation_rule = let recouperation_rate_aux=Aux{:recouperation_rate}()
        Neighbors{:endemic_presence}(Moore(2)) do data, hood, prescences, I
            any(prescences) || return prescences
            recouperation_rate = DG.get(data, recouperation_rate_aux)
            nbr_sums = sum(hood)
            map(prescences, nbr_sums, recouperation_rate) do p, n_nbrs, rr
                if p
                    true
                else
                    rand(Float32) < (n_nbrs * rr / length(hood))
                end
            end
        end
    end

    risks_rule = ExtirpationRisks{:endemic_presence}(;
        f=tanh,
        pred_pop=Grid{:pred_pop}(),
        pred_suscept=Aux{:pred_suscept}(),
        stochastic_extirpation=0.005f0/aggfactor,
    )

    aux_pred_risks_rule = ExtirpationRisks{:endemic_presence}(;
        f=tanh,
        pred_pop=Aux{:pred_pop}(),
        pred_suscept=Aux{:pred_suscept}(),
        stochastic_extirpation=0.001f0/aggfactor,
    )

    tspans = map(introductions) do intros
        # t1 = minimum(x -> x.year, intros) - 2
        1549:2018
    end
    inits = map(pred_pops, pred_carrycaps, endemic_presences) do pred_pop, pred_carrycap, endemic_presence
        (; pred_pop, pred_carrycap, endemic_presence)
    end
    pred_inits = map(pred_pops, pred_carrycaps) do pred_pop, pred_carrycap
        (; pred_pop, pred_carrycap)
    end
    endemic_inits = map(endemic_presences) do endemic_presence
        (; endemic_presence)
    end

    pred_ruleset = Ruleset(
        pred_carrycap_rule,
        introduction_rule,
        pred_spread_rule,
        pred_growth_rule;
        boundary=Remove()
    )

    endemic_ruleset = Ruleset(
        endemic_recouperation_rule,
        aux_pred_risks_rule,
        clearing_rule;
        boundary=Remove()
    )

    ruleset = Ruleset(
        DynamicGrids.rules(pred_ruleset)...,
        endemic_recouperation_rule, risks_rule, clearing_rule;
        boundary=Remove()
    )

    rules = (;
        introduction_rule,
        pred_carrycap_rule,
        pred_spread_rule,
        pred_growth_rule,
        endemic_recouperation_rule,
        risks_rule,
        aux_pred_risks_rule,
        clearing_rule,
    )

    outputs_kw = map(island_names, ns_extinct, tspans, auxs, pred_pops_aux) do island, n_extinct, tspan, aux1, pred_pop
        aux = (;
            introductions=getproperty(introductions, island),
            pred_suscept=getproperty(pred_suscepts, island),
            recouperation_rate=getproperty(recouperation_rates, island),
            habitat_requirement=getproperty(habitat_requirements, island),
            dem=getproperty(stencil_dems, island),
            pred_pop,
            landcover=aux1.lc
        )
        (; aux, mask=getproperty(stencil_masks, island), replicates, tspan, n_extinct)
    end
    outputs = map(inits, outputs_kw) do init, kw
        ResultOutput(init; kw...)
    end
    max_outputs = map(inits, outputs_kw) do init, kw
        TransformedOutput(init; kw...) do f
            presentces = Base.reduce(parent(parent(f.endemic_presence)); init=zero(eltype(f.endemic_presence))) do acc, xs
                map(|, acc, xs)
            end
            mean(f.pred_pop)
        end
    end
    pred_outputs = map(pred_inits, outputs_kw) do init, kw
        if isnothing(replicates)
            trans_output = TransformedOutput(init; kw...) do f
                Array(f.pred_pop)
                # Base.reduce(parent(parent(f.endemic_presences)); init=zero(eltype(f.endemic_presences))) do acc, xs
                #     map(|, acc, xs)
                # end
            end
        else
            ResultOutput(init; kw...)
        end
    end

    islands = map(inits, endemic_inits, pred_inits, outputs, max_outputs, pred_outputs, outputs_kw, auxs) do init, endemic_init, pred_init, output, max_output, pred_output, output_kw, aux
        (; init, endemic_init, pred_init, output, max_output, pred_output, output_kw, mask=aux.mask)
    end

    return (; ruleset, rules, pred_ruleset, endemic_ruleset, islands)
end
