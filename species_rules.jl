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

struct InteractiveCarryCap{R,W,CC,CS,I<:NamedTuple} <: Dispersal.GrowthRule{R,W}
    carrycap::CC
    carrycap_scaling::CS
    inputs::I
end
function InteractiveCarryCap{R,W}(; carrycap, carrycap_scaling, inputs=(;),) where {R,W}
    InteractiveCarryCap{R,W}(carrycap, carrycap_scaling, inputs)
end

function DynamicGrids.applyrule(
    data, rule::InteractiveCarryCap, 
    populations::NamedVector{Keys}, I,
) where Keys
    # Somehow the closure breaks without this
    g = let data=data, I=I
        val -> (data, val, I)
    end
    local_inputs = map(g, rule.inputs)
    relative_pop = NamedTuple(populations ./ rule.carrycap)
    # combine populations
    params = merge(relative_pop, local_inputs)
    
    scaling = map(rule.carrycap_scaling) do val_f
        f = DynamicGrids._unwrap(val_f)
        oneunit(eltype(populations)) + f(params)
    end |> NamedVector
    # # Carrycap cant be zero
    new_carrycaps = max.(oneunit(eltype(rule.carrycap)) .* 1f-10, rule.carrycap .* scaling)
    # any(isinf, new_carrycaps) && error("Inf carrycap found: $new_carrycaps from relative_pop $relative_pop params $params populations $populations, and scaling $scaling")
    # any(isnan, local_inputs) && error("""
    #     NaN carrycap found: $new_carrycaps 
    #         from relative_pop $relative_pop 
    #         params $params 
    #         populations $populations 
    #         and scaling $scaling
    # """)
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
    pred_suscept = get(data, rule.pred_suscept)[1]
    map(endemic_presence, pred_suscept) do present, ps
        if present
            # rand() < hp * hs & rand() < sum(map(*, ps, pred_pop))
            rand(Float32) > rule.f(sum(map(*, ps, pred_pop))) || rand(Float32) > rule.stochastic_extirpation
            # ... etc
        else
            false
        end
    end
end

    # # scaling = scale_carrycap(populations, rule.carrycap, rule.interactions, local_suitabilities)

    # newpops = map(populations, growth_rates, final_carrycap) do N, rt, cc
    #     if rt > zero(rt)
    #         if cc <= zero(cc)
    #             zero(N)
    #         else
    #             (N * cc) / (N + (cc - N) * exp(-rt))
    #         end
    #     else
    #         N * exp(rt)
    #     end
    # end
    # any(isnan, newpops) && error("NaN population found")
    # out = max.(zero(eltype(populations)), newpops)
    # # any(map(>, out, final_carrycap)) && @show out rule.carrycap final_carrycap
    # return out
# end

function def_syms(
    pred_df, introductions_df, aggfactor, dems, masks, slope_stacks, island_extinct_tables, auxs;
    replicates=nothing, pred_pops_aux=map(_ -> nothing, dems),
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
        cat =        p -> 0.0f0,#1.0p.black_rat + 0.2p.norway_rat + 1.0p.mouse + 10p.urban + 2p.cleared,
        black_rat  = p -> 0.0f0,#-0.2p.cat - 0.1p.norway_rat - 0.1p.mouse + 0.5p.native + 0.3p.abandoned + 1p.urban,
        norway_rat = p -> 0.0f0,#-0.1p.cat - 0.1p.black_rat - 0.1p.mouse + 1.5p.urban - 0.2p.native,
        mouse =      p -> 0.0f0,#-0.3p.cat - 0.2p.black_rat - 0.2p.norway_rat + 0.8p.cleared + 1.5p.urban,
        pig =        p -> 0.0f0,#0.5p.native + 0.4p.abandoned - 1p.urban - 0.3p.cleared,
        snake =      p -> 0.0f0,#-0.2p.cat + 0.2p.black_rat + 0.3p.mouse - 0.5p.urban + 0.3p.native,
        macaque =    p -> 0.0f0,#4p.abandoned + 2p.forestry + 1.5p.native - 1.0p.urban - 0.9p.cleared
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

    aux_keys = keys(first(auxs))
    pred_inputs = NamedTuple{aux_keys}(map(Aux, aux_keys))

    # How much these species are supported outside of this system
    # populations = NV(cat=0.01, black_rat=25.0, norway_rat=10.0, pig=0.3)

    # suitabilities = (human_dependency => human_intensity, forest_preference => forest_density)
    # scale_carrycap(populations, carrycap, interactions, suitabilities)

    ag_masks = map(masks) do mask
        Rasters.aggregate(Rasters.Center(), mask, aggfactor)
    end
    ag_slope = map(slope_stacks) do s
        Rasters.aggregate(maximum, replace_missing(s.slope, 0), aggfactor)
    end
    ag_auxs = map(auxs) do island_aux
        ag_aux = map(island_aux) do A
            ag_A1 = DimensionalData.modify(BitArray, Rasters.aggregate(Center(), A, (X(aggfactor), Y(aggfactor))))
            # Extend a century earlier
            ag_A = Rasters.extend(ag_A1; to=(Ti(Sampled(1500:1:2020; sampling=Intervals(Start())))))
            broadcast_dims!(identity, view(ag_A, Ti=1500..1600), view(ag_A1, Ti=At(1600)))
            return ag_A
        end
    end

    ag_roughness = map(dems) do dem
        # Clip roughness at a maximum of 500
        rgh = min.(roughness(replace_missing(dem, 0)), 500)
        Float32.(replace_missing(Rasters.aggregate(maximum, rgh ./ 500, aggfactor), NaN))
    end
    ag_dems = map(dems) do dem
        Float32.(replace_missing(Rasters.aggregate(Rasters.Center(), dem, aggfactor), NaN))
    end
    ag_local_cost = map(ag_dems) do dem
        window = Window{2}()
    end

    kernel = DispersalKernel(
        stencil=moore,
        formulation=ExponentialKernel(Param(1.0f0, bounds=(0.0f0, 2.0f0))),
        cellsize=1.0f0,
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
        1.0f0
    end

    introductions = map(island_names) do key
        island_df = filter(r -> r.island == string(key), introductions_df)
        inits = map(pred_indices, one_individual) do i, oi
            x = zeros(Float16, length(pred_keys))
            x[i] = Float16(oi * 100 * aggfactor)
            PredNV(x)
        end

        map(eachrow(island_df)) do r
            init = inits[Symbol(r.species)]
            (; year=r.year, geometry=(X=r.lon, Y=r.lat), init)
        end
    end

    pred_rmax = Float32.(PredNV(pred_df.rmax))
    pred_max_density = Float32.(PredNV(pred_df.max_density))
    pred_pops = map(ag_masks) do m
        map(_ -> map(_ -> 0.0f0, pred_rmax), m)
    end
    pred_carrycaps = map(ag_masks) do m
        map(_ -> carrycap, m)
    end

    pred_carrycap_rule = InteractiveCarryCap{:pred_pop,:pred_carrycap}(;
        carrycap,
        carrycap_scaling=map(Val, pred_funcs),
        inputs=pred_inputs,
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
        cat_sus = ExtinctNV(LinRange(Float32(0.0), Float32(0.02), n_extinct))
        black_rat_sus = ExtinctNV(LinRange(Float32(0.0), Float32(0.01), n_extinct))
        norway_rat_sus = ExtinctNV(LinRange(Float32(0.0), Float32(0.01), n_extinct))
        mouse_sus = ExtinctNV(LinRange(Float32(0.0), Float32(0.001), n_extinct))
        pig_sus = ExtinctNV(LinRange(Float32(0.0), Float32(0.005), n_extinct))
        snake_sus = ExtinctNV(LinRange(Float32(0.0), Float32(0.005), n_extinct))
        macaque_sus = ExtinctNV(LinRange(Float32(0.0), Float32(0.003), n_extinct))
        pred_suscept = map(cat_sus, black_rat_sus, norway_rat_sus, mouse_sus, pig_sus, snake_sus, macaque_sus) do c, br, nr, m, p, s, ma
            PredNV((c, br, nr, m, p, s, ma)) ./ aggscale
        end
    end

    # Every species is everywhere initially, in this dumb model
    endemic_presences = map(ag_masks, hunting_suscepts) do mask, hunting_suscept
        map(mask) do m
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

    clearing_rule = let native=Aux{:native}()
        Cell{:endemic_presence}() do data, presences, I
            uc = get(data, native, I)
            presences .& uc
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
    pred_spread_rule = let demaux=Aux{:dem}(), aggfactor=aggfactor, pred_kernels=pred_kernels, carrycap=carrycap
        SetNeighbors{Tuple{:pred_pop,:pred_carrycap}}(pred_kernels[1]) do data, hood, (Ns, _), I
            if any(isnan, Ns) 
                data[:pred_pop][I...] = zero(Ns)
                return nothing
            end
            Ns == zero(Ns) && return nothing
            dem = get(data, demaux)
            reps = DynamicGrids.replicates(data)
            carrycap_nbrs = if isnothing(reps)
                DG.neighbors(DG.grids(data).pred_carrycap, I)
            else
                DG.neighbors(DG.grids(data).pred_carrycap, (I..., reps))
            end
            elev_nbrs = DG.neighbors(dem, I)
            @inbounds elev_center = dem[I...]

            # spreadability = DG.neighbors(DG.mask(data), I)
            sum = zero(Ns) # TODO needs cell size here
            cellsize = (30 * aggfactor)

            # Randomise hood starting position to avoid
            # directional artifacts in output
            start = rand(0:length(hood)-1)
            @inbounds for ix in eachindex(hood)
                # Rotate indices in relation to starting point
                i = start + ix
                if i > length(hood)
                    i = i - length(hood)
                end
                sp = carrycap_nbrs[i] ./ carrycap
                any(map(isnan, sp)) && continue
                Ih = DG.indices(hood, I)[i]
                ks = getindex.(DynamicGrids.kernel.(pred_kernels), i)
                d = DG.distances(hood)[i]
                e = elev_nbrs[i]
                # any(isnan, N) && error("N is NaN")
                # isnan(k) && error("k is NaN")
                # isnan(sp) && error("sp is NaN")
                # slope = (e - elev_center) / (cellsize * d)
                # Imhof Tobler equation
                # speed = (6ℯ.^(slopefactor .* (slope .+ distfactor)))

                # slope_factor = 1 - min(slope, 1) # TODO needs cell size here
                slope_factor = 1.0000000001f0 - min(abs(e - elev_center) / cellsize, 1)
                # @show N  k  sp slope_factor one_individual
                # propagules = N .* k .* sp .* slope_factor .+ ((rand(typeof(N)) .- 0.5) .* one_individual)
                # propagules = Ns * sp * slope_factor .* kr .* rand(typeof(Ns))
                propagules = trunc.(Float32.(Ns .* sp .* slope_factor .* rand(typeof(Ns)) .^ 3 .* ks .* 40))
                # any(isnan, propagules) && error("NaNs found Ns: $Ns sp: $sp slope_factor: $slope_factor ks: $ks")
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
            # @inbounds any(isnan, data[:pred_pop][I...]) && error("NaN in output")
            @inbounds sub!(data[:pred_pop], sum, I...)
            return nothing
        end
    end # let

    # This is not really allee extinction, just assuring
    # at least one individual exists.
    # pred_allee = AlleeExtinction{:pred_pop}(;
    #     minfounders=one_individual,
    # )

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
        stochastic_extirpation=0.001f0/aggfactor, 
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

    pred_ruleset = Ruleset(
        introduction_rule,
        pred_carrycap_rule, 
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
        Chain(endemic_recouperation_rule, risks_rule, clearing_rule);
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

    outputs_kw = map(island_names, ns_extinct, tspans, ag_auxs, pred_pops_aux) do island, n_extinct, tspan, ag_aux, pred_pop
        aux = (;
            spreadability=getproperty(spreadability, island),
            introductions=getproperty(introductions, island),
            pred_suscept=[getproperty(pred_suscepts, island)],
            recouperation_rate=getproperty(recouperation_rates, island),
            habitat_requirement=getproperty(habitat_requirements, island),
            dem=getproperty(stencil_dems, island),
            pred_pop,
            ag_aux...
        )
        (; aux, mask=getproperty(stencil_masks, island), replicates, tspan, n_extinct)
    end
    outputs = map(inits, outputs_kw) do init, kw
        ResultOutput(init; kw...)
        # TransformedOutput(init; kw...) do f
        #     presentces = Base.reduce(parent(parent(f.endemic_presences)); init=zero(eltype(f.endemic_presences))) do acc, xs
        #         map(|, acc, xs)
        #     end
        #     Base.maximum(f.pred_pop)
        # end
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
    # Sum as we go
    # ncells = sum(DG.mask(output))
    # last_obs = map(island_names) do island
    #     ExtinctNV(map(first, getproperty(island_extinct_tables, island).extinct))
    # end

    islands = map(inits, pred_inits, outputs, pred_outputs, outputs_kw, ag_masks) do init, pred_init, output, pred_output, output_kw, ag_mask 
        (; init, pred_init, output, pred_output, output_kw, ag_mask)
    end


    return (; ruleset, rules, pred_ruleset, endemic_ruleset, islands)
end
