pn(x, d) = Distributions.pdf(x, d) ./ Distributions.pdf(x, 0)

# The amount of influence neighbors have
b = (; bounds=(-2.0, 2.0))
transition_probability = let empty = (native=0.0, cleared=0.0, abandoned=0.0, urban=0.0, forestry=0.0, water=0.0,),
                  k = Exponential(P(2.0; bounds=(0.00001, 5.0), label="kernel")),
                  cn = P(0.1; b..., label="cleared from native"),
                  cc = P(1.0; b..., label="cleared from cleared"),
                  ca = P(0.2; b..., label="cleared from abandoned"),
                  aa = P(0.8; b..., label="abandoned from abandoned"),
                  uc = P(0.2; b..., label="urban from cleared"),
                  uu = P(1.0; b..., label="urban from urban")
    # ff = Exponential(P(1.7; b..., label="forestry from forestry"))
    (
        native = empty,
        cleared = (
            native=d -> cn * pn(k, d), 
            cleared=d -> cc * pn(k, d), 
            abandoned=d -> ca * pn(k, d), 
            urban=0.0,
            forestry=0.0, 
            water=0.0,
        ),
        abandoned = (
            native=0.0,
            cleared=0.0,
            abandoned=d -> aa * pn(k, d),
            urban=0.0,
            forestry=0.0,
            water=0.0,
        ),
        urban = (
            native=0.0, 
            cleared=d -> uc * pn(k, d), 
            abandoned=0.0, 
            urban=d -> uu * pn(k, d), 
            forestry=0.0, 
            water=0.0,
        ),
        forestry = empty,
        water = empty,
    )
end

# The logic of sequential category change - can a transition happen at all
# Human Population and species introduction events
eventrule = let states=states
    # D=dims(masks.mus),
    SetGrid{:landcover,:landcover}() do data, l1, l2
        history = DynamicGrids.aux(data).history
        events = DynamicGrids.aux(data).events
        current_year = currenttime(data)
        if hasselection(history, Ti(At(current_year)))
            foreach(eachindex(l1), l1, view(history, Ti(At(current_year)))) do I, state, hist
                # Fill water
                # if hist.water && count(hist) == 1
                    # l1[I] = states.water 
                # end
                if count(hist) == 1
                    l1[I] = findfirst(hist)
                end
            end
        end
    end
end
# for event in events
#     if event.year == current_year
#         p = event.geometry
#         I = DimensionalData.dims2indices(D, (X(Contains(p.X)), Y(Contains(p.Y))))
#         Iu = map(i -> i-1:i+1, I)
#         Ic = map(i -> i-10:i+10, I)
#         l1[Iu...] .= states.urban
#         l1[Ic...] .= states.cleared
#         l2[Iu...] .= states.urban
#         l2[Ic...] .= states.cleared
#         # l2 .= states.cleared
#     end
# end
        
# pressure = NamedVector(
#     native=P(1.0; b...),
#     cleared=P(1.5; b...),
#     abandoned=P(1.0; b...),
#     urban=P(1.0; b...),
#     forestry=P(1.0; b...),
# )
pressure = let preds=lc_targets.mus, ngridcells=size(sum(masks.mus))
    leverage=P(3.0; bounds=(1.0, 10.0))
    # cleared=P(1.5; b...),
    # urban=P(1.4; b...)
    abandoned = 0.0
    forestry = Inf
    native = 0.0 # Never regrows
    water = Inf # Put it in as soon as its on the next map
    (data, rule) -> begin
        predicted = preds[At(currenttime(data))]
        stategrid = data[:landcover]
        hist = view(DynamicGrids.aux(data).history, Ti(Contains(currenttime(data))))
        allowed_cleared = ThreadsX.count(x -> x.cleared, hist; init=0)
        allowed_urban = ThreadsX.count(x -> x.urban, hist; init=0)
        ncleared = ThreadsX.sum(==(rule.states.cleared), parent(stategrid); init=0)
        nurban = ThreadsX.sum(==(rule.states.urban), parent(stategrid); init=0)
        urban = LandscapeChange.calc_pressure(leverage, nurban, predicted.urban, allowed_urban)
        cleared = LandscapeChange.calc_pressure(leverage, ncleared, predicted.cleared, allowed_cleared)
        # @show cleared urban ncleared nurban ncleared allowed_urban allowed_cleared
        NamedVector(; native, cleared, abandoned, urban, forestry, water)
    end
end
# A = map(xy -> calc_pressure(1, 0.1, xy..., 10000), DimPoints((n=0:100, p=0:100)))
# Makie.heatmap(A)
# v = 1.0 + (-log(rand()^1.0))

n, c, a, u, f, w = states
precursors = (
    native =    SA[n, n, n, n],
    cleared =   SA[n, c, a, f],
    abandoned = SA[c, a, a, f],
    urban =     SA[n, c, a, f],
    forestry =  SA[n, c, a, f],
    water =     SA[w, w, w, w],
)

staterule = BottomUp{:landcover}(;
    stencil=Moore(2),
    states ,
    inertia=P(1.0),
    logic=LandscapeChange.all_transitions(transitions),
    transitions=transition_probability,
    suitability=map(_ -> 1, states), #Aux{:suitability}(),
    history=Aux{:history}(),
    fixed=false,
    perturbation=P(2.0; bounds=(0.0, 10.0), label="perturbation"),
)

degradationrule = let states=states, 
                      degradation_curve=Param(1.0; bounds=(0.0, 4.0)), 
                      degradation_rate=Param(0.01; bounds=(0.0, 0.01))
    kernl = Kernel(Moore{4}()) do d
        exp(-d / degradation_curve)# * degredation_rate
    end
    Neighbors{Tuple{:native_fraction,:landcover},:native_fraction}(; stencil=kernl) do data, hood, (native_fraction, landcover), I
        DynamicGrids.ismasked(data, I...) && return oneunit(native_fraction)
        target_native_fraction = get(data, Aux{:target_native_fraction}(), I)
        new_native_fraction = if landcover == states.native
            degradation = zero(first(hood))
            @simd for i in 1:length(hood)
                @inbounds degradation += (1 - hood[i]) * kernel(hood)[i] * degradation_rate
            end
            max(native_fraction / (1 + degradation), zero(native_fraction))
        else
            zero(native_fraction)
        end
        # We don't want to overshoot the real native fraction in real final run
        return max(target_native_fraction, new_native_fraction)
    end
end

ruleset = Ruleset(staterule, eventrule; proc=CPUGPU());
