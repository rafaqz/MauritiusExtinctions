using DynamicGrids, Dispersal, LandscapeChange, ModelParameters
using StaticArrays
using Rasters
using ThreadsX
using Distributions
using Setfield
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

@inline Base.@assume_effects :foldable function DynamicGrids.applyrule(
    data, rule::InteractiveCarryCap,
    populations::NamedVector{Keys}, I,
) where Keys
    local_inputs = NamedTuple(get(data, rule.inputs, I))
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

struct ExtirpationRisks{R,W,F,T,PR,PS,PP,E,SE} <: DynamicGrids.CellRule{R,W}
    f::F
    traits::T
    pred_response::PR
    pred_suscept::PS
    pred_pop::PP
    pred_effect::E
    stochastic_extirpation::SE
end
function ExtirpationRisks{R,W}(; f, traits, pred_response, pred_suscept, pred_pop, pred_effect, stochastic_extirpation) where {R,W}
    ExtirpationRisks{R,W}(f, traits, pred_response, pred_suscept, pred_pop, pred_effect, stochastic_extirpation)
end

@inline function DynamicGrids.applyrule(data, rule::ExtirpationRisks, endemic_presence, I)
    # If they are all absent do nothing
    any(endemic_presence) || return endemic_presence

    # Get the effect of predators on each endemic
    pred_effect = if isnothing(rule.pred_effect)
        pred_pop = get(data, rule.pred_pop, I)
        pred_suscept = get(data, rule.pred_suscept)
        map(rule.f, predator_effect(pred_pop, pred_suscept))
    else
        # We have already precalculated the effect
        get(data, rule.pred_effect, I)
    end

    return map(endemic_presence, pred_effect) do present, effect
        if present 
            rand(typeof(effect)) > effect && rand() > rule.stochastic_extirpation
            # typeof(effect)(0.5) > effect
        else
            false
        end
    end
end

Base.@assume_effects :foldable function predator_effect(pred_pop, pred_suscept)
    mapreduce(+, pred_suscept, pred_pop) do ps, pp
         map(xs -> map(Float32 ∘ *, xs, pp), ps)
    end
end

# function predator_suceptibility(mass_response, pred_response, traits)
#     pred_suscept = mapreduce(+, pred_response, traits) do pr, t
#         map(mass_response, pr) do m, p
#             t .* p .* m
#         end
#     end ./ (32 * 8^2)
# end

function predator_suceptibility(pred_response, traits)
    mapreduce(+, pred_response, traits) do pr, t
        map(NamedVector(pr)) do p
            NamedVector(t) .* p
        end |> NamedVector
    end ./ (32 * 8^2)
end

@inline function DynamicGrids.modifyrule(rule::ExtirpationRisks, data::AbstractSimData)
    if isnothing(rule.pred_effect) 
        pred_response = get(data, rule.pred_response)
        # mass_response = get(data, rule.mass_response)
        traits = get(data, rule.traits)
        @set! rule.pred_suscept = predator_suceptibility(pred_response, traits)
    end
    @set! rule.traits = nothing # Simplify for GPU argument size
    @set rule.pred_response = nothing # Simplify for GPU argument size
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
    @show typeof(lc_ag)
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
    pred_df, introductions_df, island_endemic_tables, auxs, aggfactor;
    replicates=nothing, pred_pops_aux,
    pred_keys=(:cat, :black_rat, :norway_rat, :mouse, :pig, :wolf_snake, :macaque)
)
    aggscale = aggfactor^2
    window = Window{2}()
    moore = Moore{3}()
    ns_endemic = map(nrow, island_endemic_tables)
    endemic_keys = map(island_endemic_tables) do endemic_table
        Tuple(Symbol.(replace.(endemic_table.Species, Ref(' ' => '_'))))
    end
    EndemicNVs = map(endemic_keys) do ek
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

    # These are taken from the literature in contexts where it seems also applicable to the Mascarenes
    carrycap = Float32.(NV(;
        cat =        0.02,
        black_rat =  30.0,
        norway_rat = 15.0, # This one is more of a guess
        mouse =      52.0,
        pig =        0.3,
        wolf_snake = 20.0, # As is this one
        macaque =    1.0,
    ) .* aggscale)[pred_keys]

    # These are parametries from reading qualitative literature and quantitative literature without
    # enough context to really use the number.
    # Maybe we should fit them to some observations in specific points at specific times.
    # Refs: Smucker et al 2000 - Hawaiii cats, rats, mice
    pred_funcs = (;
        cat =        p -> 1.0f0p.black_rat + 0.3f0p.norway_rat + 1.0f0p.mouse + 10f0p.urban + 2f0p.cleared,
        black_rat  = p -> -0.2f0p.cat - 0.1f0p.norway_rat - 0.1f0p.mouse + 0.5f0p.native + 0.3f0p.abandoned + 1p.urban,
        norway_rat = p -> -0.1f0p.cat - 0.1f0p.black_rat - 0.1f0p.mouse + 1.5f0p.urban - 0.2f0p.native,
        mouse =      p -> -0.3f0p.cat - 0.2f0p.black_rat - 0.2f0p.norway_rat + 0.8f0p.cleared + 1.5f0p.urban,
        pig =        p -> 0.5f0p.native + 0.4f0p.abandoned - 2f0p.urban - 1.0f0p.cleared,
        wolf_snake =      p -> -0.2f0p.cat + 0.2f0p.black_rat + 0.3f0p.mouse - 0.5f0p.urban + 0.3f0p.native,
        macaque =    p -> 1.0f0p.abandoned + 0.7f0p.forestry + 0.4f0p.native - 1.0f0p.urban - 0.8f0p.cleared
    )[pred_keys]

    # These need to somewhat balance low growth rates. They are almost totally made up.
    # The units are in pixels - it needs fixing to the aggregation size.
    spread_rate = NV(;
        cat =         30.0f0,
        black_rat =   1.0f0,
        norway_rat =  1.0f0,
        mouse =       0.5f0,
        pig =         15.0f0,
        wolf_snake =  1.0f0,
        macaque =     5.0f0,
    )

    # These dummy values will be optimised. Its likely too many for that but we
    # can reduce the numver of predators, removing mouse, macaque etc
    pred_response_raw = (;
        ismammal = (;
            cat =         0.1,
            black_rat =   0.02,
            norway_rat =  0.02,
            mouse =       0.00,
            pig =         0.04,
            wolf_snake =  0.02,
            macaque =     0.02,
        ),
        isbird = (;
            cat =         0.1,
            black_rat =   0.03,
            norway_rat =  0.03,
            mouse =       0.02,
            pig =         0.01,
            wolf_snake =  0.03,
            macaque =     0.01,
        ),
        isreptile = (;
            cat =         0.2,
            black_rat =   0.02,
            norway_rat =  0.02,
            mouse =       0.01,
            pig =         0.02,
            wolf_snake =  0.02,
            macaque =     0.02,
        ),
        isgroundnesting = (;
            cat =         0.2,
            black_rat =   0.01,
            norway_rat =  0.02,
            mouse =       0.01,
            pig =         0.05,
            wolf_snake =  0.0,
            macaque =     0.0,
        ),
        flightlessness = (;
            cat =         0.2,
            black_rat =   0.01,
            norway_rat =  0.02,
            mouse =       0.01,
            pig =         0.05,
            wolf_snake =  0.01,
            macaque =     0.0,
        ),
    )

    response_keys = NamedTuple{keys(pred_response_raw)}(keys(pred_response_raw))
    pred_response = map(response_keys, pred_response_raw) do k, ps
        map(NamedTuple(ps[pred_keys]), NamedTuple{pred_keys}(pred_keys)) do val, label
            Param(val; label=Symbol(k, :_, label))
        end
    end

    island_endemic_traits = map(island_endemic_tables, ns_endemic, EndemicNVs) do table, n_endemic, EndemicNV
        ismammal = EndemicNV(table.Group .== "mammal")
        isbird = EndemicNV(table.Group .== "bird")
        isreptile = EndemicNV(table.Group .== "reptile")
        isgroundnesting = EndemicNV(table.Ground_nesting)
        flightlessness = EndemicNV(Float32.(1 .- table.Flight_capacity))
        (; ismammal, isbird, isreptile, isgroundnesting, flightlessness)
    end
    gecko_mass = 8 # estimated mean of multiple species
    skink_mass = 3 # estimated mean of multiple species
    mouse_mass = 16.25
    wolf_snake_pre_mass = round(0.48gecko_mass + 0.30mouse_mass + 0.22skink_mass)

    mean_prey_mass = (;
        cat =         (41.0, 51.0), # Pearre and Maaas 1998
        black_rat =   (10.0, 10.0), # made up
        norway_rat =  (8.0f0, 8.0), # made up
        mouse =       (3.0f0, 5.0), # made up
        pig =         (100.0, 100.0), # made up
        wolf_snake =  (9.0f0, 7.0), # Estimated from Fritts 1993
        macaque =     (100.0, 100.0), # made up
    )

    island_mass_response = map(island_endemic_tables, EndemicNVs) do table, EndemicNV
        endemic_mass = EndemicNV(table.Mass)
        map(mean_prey_mass) do (mean, std)
            dist = Distributions.Normal(mean, 2std) # Doubled because of the skew
            scalar = 1 / pdf(dist, mean)
            map(endemic_mass) do pm
                pdf(dist, pm) * scalar
            end
        end
    end

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
    pred_init_nvs = map(pred_indices) do i
        x = zeros(Float16, length(pred_keys))
        x[i] = 50 * aggfactor
        PredNV(x)
    end

    introductions = map(island_names) do key
        island_df = filter(r -> r.island == string(key), introductions_df)
        map(eachrow(island_df)) do r
            init = pred_init_nvs[Symbol(r.species)]
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

    habitat_requirements = map(EndemicNVs) do EndemicNV
        Float32.(rand(EndemicNV))
    end

    hunting_suscepts = map(EndemicNVs) do EndemicNV
        Float32.(rand(EndemicNV))
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
            presences .& (lc.native > 0.0f0)
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
    pred_spread_rule = let demaux=Aux{:dem}(), aggfactor=aggfactor, pred_kernels=pred_kernels, carrycap=carrycap, slope_scalar=5.0 # Param(5.0, bounds=(1.0, 20.0))
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

    recouperation_rates = map(EndemicNVs) do EndemicNV
        Float32.(rand(EndemicNV) / 10)
    end
    endemic_recouperation_rule = let recouperation_rate_aux=Aux{:recouperation_rate}()
        Neighbors{:endemic_presence}(Moore(1)) do data, hood, prescences, I
            # any(prescences) || return prescences
            recouperation_rate = DG.get(data, recouperation_rate_aux)
            nbr_sums = foldl(hood; init=Base.reinterpret.(UInt8, zero(first(hood)))) do x, y
                Base.reinterpret(UInt8, x) + Base.reinterpret(UInt8, y)
            end
            map(prescences, nbr_sums, recouperation_rate) do p, n_nbrs, rr
                if p
                    true
                elseif n_nbrs > 0
                    rand(Float32) < (n_nbrs * rr / length(hood))
                else
                    false
                end
            end
        end
    end

    risks_rule = ExtirpationRisks{:endemic_presence}(;
        f=tanh,
        traits=Aux{:endemic_traits}(),
        pred_response,
        pred_suscept=nothing,
        pred_effect=nothing,
        pred_pop=Grid{:pred_pop}(),
        stochastic_extirpation=0.005f0/aggfactor,
    )

    aux_pred_risks_rule = ExtirpationRisks{:endemic_presence}(;
        f=tanh,
        traits=Aux{:endemic_traits}(),
        pred_response,
        pred_suscept=nothing,
        pred_effect=Aux{:pred_effect}(),
        pred_pop=nothing,
        stochastic_extirpation=0.005f0/aggfactor,
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
        Chain(endemic_recouperation_rule, aux_pred_risks_rule, clearing_rule);
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

    pred_effects = map(pred_pops_aux, island_mass_response, island_endemic_traits) do pred_pop, mass_response, traits
        pred_suscept = predator_suceptibility(pred_response, traits)
        isnothing(pred_pop) ? nothing : generate_predator_effect(risks_rule.f, pred_pop, pred_suscept)
    end

    outputs_kw = map(island_names, ns_endemic, tspans, auxs, pred_pops_aux, island_endemic_traits, pred_effects) do island, n_endemic, tspan, aux1, pred_pop, endemic_traits, pred_effect
        aux = (;
            introductions=getproperty(introductions, island),
            recouperation_rate=getproperty(recouperation_rates, island),
            habitat_requirement=getproperty(habitat_requirements, island),
            dem=getproperty(stencil_dems, island),
            pred_pop,
            endemic_traits,
            pred_effect,
            landcover=aux1.lc
        )
        (; aux, mask=getproperty(stencil_masks, island), replicates, tspan, n_endemic)
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
            end
        else
            ResultOutput(init; kw...)
        end
    end
    endemic_outputs = map(endemic_inits, outputs_kw) do init, kw
        ResultOutput(init; kw...)
    end
    outputs = map(inits, outputs_kw) do init, kw
        ResultOutput(init; kw...)
    end

    island_keys = NamedTuple{keys(inits)}(keys(inits))
    islands = map(
        island_keys, inits, endemic_inits, pred_inits, outputs, max_outputs, endemic_outputs, pred_outputs, outputs_kw, auxs, island_mass_response
    ) do key, init, endemic_init, pred_init, output, max_output, endemic_output, pred_output, output_kw, aux, mass_response
        (; key, init, endemic_init, pred_init, output, max_output, endemic_output, pred_output, output_kw, aux, mass_response)
    end

    return (; ruleset, rules, pred_ruleset, endemic_ruleset, islands, pred_response)
end

function gpu_cleanup(A)
    x, y, ti = lookup(A, (X, Y, Ti))
    Rasters.set(A,
        Ti => Sampled(first(ti) - 0.5:last(ti) - 0.5; sampling=Intervals(Start())),
        X => LinRange(first(x), last(x), length(x)),
        Y => LinRange(first(y), last(y), length(y)),
    )
end

function generate_predator_effect!(f, x, pred_pop, pred_suscept)
    ThreadsX.map!(x, pred_pop) do pop
        map(f, predator_effect(pop, pred_suscept))
    end
end

function generate_predator_effect(f, pred_pop::Union{AbstractArray{<:Any,2},AbstractArray{<:Any,3}}, pred_suscept)
    xs = ThreadsX.map(pred_pop) do pop
        map(f, predator_effect(pop, pred_suscept))
    end
    rebuild(pred_pop, xs)
end

function predict_extinctions(endemic_ruleset::Ruleset, islands, pred_response)
    map(islands) do island
        _predict_extinctions(endemic_ruleset, island, pred_response)
    end
end
function _predict_extinctions(endemic_ruleset::Ruleset, island, pred_response)
    @show island.key
    kw = island.output_kw
    aux = kw.aux
    pred_suscept = predator_suceptibility(pred_response, aux.endemic_traits)
    (; pred_pop, pred_effect) = aux
    generate_predator_effect!(tanh, pred_effect, pred_pop, pred_suscept)
    output = TransformedOutput(island.endemic_init; kw...) do f
        # Take the mean over the replicates dimension
        mean(sum, eachslice(f.endemic_presence; dims=3))
    end
    sim!(output, endemic_ruleset; printframe=false, proc=CPUGPU())
    return output
    # Return NamedVector of years to extinction
    # return sum(output) .+ first(island.output_kw.tspan) .- 1
end
