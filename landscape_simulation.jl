using LandscapeChange
using StatsPlots
using Unitful
using StaticArrays
# using ReverseStackTraces
using GLMakie
using Statistics, StatsBase, GLM
using ThreadsX
using ModelParameters
using DimensionalData, GLMakie
using NCDatasets, Rasters

# From woods and forests in Mauritius, p 40:40
# 1880: 70,000 acres out of 300,000 remain
# 35,000 were native (but "dilatipated and ruined?")
# Also mentions that invasives replace natives
include("travel_cost.jl")
include("tabular_data.jl")
include("landcover_logic.jl")
include("map_file_list_2.jl")

states = NamedVector(lc_categories)

const P = RealParam 
const NV = NamedVector

# to category from category
logic = NV(
    native    = NV(native=true,  cleared=false, abandoned=false, urban=false, forestry=false,  water=false),
    cleared   = NV(native=true,  cleared=true,  abandoned=true,  urban=false, forestry=false,  water=false),
    abandoned = NV(native=false, cleared=true,  abandoned=true,  urban=false, forestry=false,  water=false),
    urban     = NV(native=true,  cleared=true,  abandoned=true,  urban=true,  forestry=false,  water=false),
    forestry  = NV(native=true,  cleared=true,  abandoned=true,  urban=false, forestry=true,   water=false),
    water     = NV(native=true,  cleared=true,  abandoned=true,  urban=true,  forestry=false,  water=true),
)

slices = make_raster_slices(masks, lc_categories; category_names);
nv_rast = Rasters.combine(namedvector_raster.(slices.mus.timeline))
striped_raw = stripe_raster(nv_rast, states)
# Rasters.rplot(striped_raw[Ti=1750..2050]; colorrange=(1, 6))
compiled = compile_timeline(logic, nv_rast; assume_continuity=true)
history = compiled.timeline
striped_compiled = stripe_raster(compiled.timeline, states)
# Rasters.rplot(striped_compiled[Ti=1750..2050]; colorrange=(1, 6))

# pixel_timeline1 = compiled.timeline[Y=Near(-20.363), X=Near(57.5015)]
# pixel_timeline = compiled.timeline[Y=Near(-20.363), X=Near(57.5015)]
# indirect=indirect_logic(logic) 
# reversed=reverse_transitions(logic) 
# reversed_indirect=reverse_transitions(indirect)
# remove_intermediate_uncertainty!(pixel_timeline, logic, reversed, indirect, reversed_indirect)
# remove_intermediate_uncertainty!(view(pixel_timeline, reverse(eachindex(pixel_timeline))), reversed, logic, reversed_indirect, indirect)
# pixel_timeline == pixel_timeline1
# pixel_timeline

# a = NV(native=true, cleared=false, abandoned=false, urban=true, forestry=false, water=false)
# b = NV(native=false, cleared=true, abandoned=false, urban=false, forestry=false, water=false)
# _merge_all_possible(a, b, logic) |> pairs
# _merge_all_possible(b, a, logic) |> pairs
logic

cat_counts = let states=states
    map(slice(history[Ti=1600..2017], Ti)) do slice
        total_counts = zeros(Int, size(first(slice)))
        known_counts = zeros(Int, size(first(slice)))
        for categories in slice
            total_counts .+= categories
            if count(categories) == 1
                known_counts .+= categories
            end
        end
        vals = ntuple(length(states)) do i
            known = known_counts[i]
            total = total_counts[i]
            year = first(refdims(slice, Ti))
            pop = human_pop_timeline.mus[Near(year)]
            (; total, known, ratio=known/total, meancount=(known + total) / 2, year, pop)
        end 
        NamedVector{keys(states)}(vals)
    end 
end

high_certainty = map(category_names) do k
    vals = map(cat_counts) do val
        val[k]
    end
    filter(v -> v.ratio > 0.7, vals)
end

# Cleared land is used for urbanisation by 1992, so don't use it in the model
cleared_model = lm(@formula(meancount ~ pop^2 + pop), DataFrame(high_certainty.cleared))
urban_model = lm(@formula(meancount ~ pop^2), DataFrame(high_certainty.urban))
ti = dims(human_pop_timeline.mus, Ti)
pops = map(pop -> (; pop), human_pop_timeline.mus)
cleared_pred = DimArray(predict(cleared_model, parent(pops)), ti)
urban_pred = DimArray(predict(urban_model, parent(pops)), ti)
lc_predictions = map((cleared, urban) -> (; cleared, urban), cleared_pred, urban_pred)
# Plots.plot(first.(pops))
# Plots.plot(cleared_pred)
# Plots.scatter!(getproperty.(high_certainty.cleared, :known))
# Plots.scatter!(getproperty.(high_certainty.cleared, :meancount))
# Plots.scatter!(getproperty.(high_certainty.cleared, :total))
# Plots.plot(urban_pred)
# Plots.scatter!(getproperty.(high_certainty.urban, :known))
# Plots.scatter!(getproperty.(high_certainty.urban, :meancount))
# Plots.scatter!(getproperty.(high_certainty.urban, :total))
# Plots.plot!(human_pop_timeline.mus)

# b = (; bounds=(0.0, 2.0))
# inertia = NamedVector(
#     native=P(0.2; b...),
#     cleared=P(0.6; b...),
#     abandoned=P(0.1; b...),
#     urban=P(0.9; b...),
#     forestry=P(0.9; b...),
#     water
# )
b = (; bounds=(1.0, 2.0))
landscape_events = (
    mus = [
        (year=1638, n=50, geometry=(X=57.7228, Y=-20.3754)), # "First Dutch settlement"
        (year=1664, n=160, geometry=(X=57.7228, Y=-20.3754)),  # "Second dutch settlement"
        (year=1723, n=524, geometry=(X=57.7228, Y=-20.3754)),  # "Colony at grand port founded", "French",
        (year=1723, n=20#=?=#, geometry=(X=57.5012, Y=-20.1597)),  # "Colony at port louie began", "French",
        (year=1735, n=20#=?=#, geometry=(X=57.5012, Y=-20.1597)),  # "Colony at port louie as capital", "French",
    ],
    reu = [
        (year=1665, n=20#=?=#, geometry=(X=55.4485, Y=20.8785)),  # "Colony at port louie as capital", "French",
    ],
    rod = [
        # TODO: put real dates here
        (year=1665, n=20#=?=#, geometry=(X=55.4485, Y=20.8785)),
    ]
)

# Rasters.rplot(ustrip.(last(travel_times.mus)))
human_suitability = map(travel_times, slope_stacks, dems) do tt, ss, dem
    slope_suitability = (1 .- ss.slope) # Steep bad, flat good
    slices = map(tt) do travel_time
        travel_suitability = (1 .- travel_time ./ 50u"hr")
        reverse(replace_missing(travel_suitability .* slope_suitability, 0.0); dims=Y)
    end
    Rasters.combine(slices, Ti)
end
# travel_times.mus ./ u"hr" |> Rasters.rplot
# extrema.(collect(skipmissing.(travel_times.mus)))
# (1 .- slope_stacks.mus.slope) |> Rasters.rplot
# Rasters.aggregate(sum, human_suitability.mus, 50; skipmissingval=true) |> Rasters.rplot

distance_to_water = map(island_keys) do k
    fix_order(Raster(joinpath(distancedir, string(k), "to_water.tif")))
end

suitability = map(human_suitability, distance_to_water) do hs, dtw
    # dtw = 1 ./ (1 .+ sqrt.(replace_missing(dtw, Inf)))
    native = fill!(similar(hs), 1.0)
    cleared = hs .^ 2# .* dtw
    abandoned = 1 .- hs
    urban = hs .^ 2# .* dtw
    forestry = fill!(similar(hs), 1.0)
    map(native, cleared, abandoned, urban) do n, c, a, u
        NamedVector(native=n, cleared=c, abandoned=a, urban=u, forestry=n, water=n)
    end
end

# Plots.plot(distance_to_water.rod)
# Plots.plot(distance_to_water.reu)
# Plots.plot(human_suitability.mus[Ti=6])
# Plots.plot(human_suitability.rod[Ti=1])
# Plots.plot(human_suitability.reu[Ti=20])
# pop_density.reu |> Plots.plot

# Check that suitability makes sense
# suit = human_suitability.rod
# pop = pop_density.rod

# pop_models = map(human_suitability, pop_density) do suit, pop
#     ag_pop_density = Rasters.aggregate(sum, pop, 20; skipmissingval=true)
#     ag_pop_density = replace(ag_pop_density, NaN=>0.0)
#     resample_suit = mask!(resample(suit[Ti=1]; to=ag_pop_density); with=ag_pop_density)
#     mask!(ag_pop_density; with=resample_suit)
#     df = RasterStack((pop=replace_missing(ag_pop_density), suit=replace_missing(resample_suit)))
#     model = lm(@formula(pop ~ suit), df)
# end
# map(r2, pop_models)
# Plots.plot(dems.mus)
# Rasters.aggregate(sum, pop_density.reu, 50; skipmissingval=true) |> Makie.plot
# Rasters.aggregate(sum, human_suitability.reu, 50; skipmissingval=true) |> rplot

pn(x, d) = Distributions.pdf(x, d) ./ Distributions.pdf(x, 0)

# The amount of influence neighbors have
b = (; bounds=(-2.0, 2.0))
transitions = let empty = (native=0.0, cleared=0.0, abandoned=0.0, urban=0.0, forestry=0.0, water=0.0,),
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
map(==(keys(first(transitions))) ∘ keys, transitions)

# The logic of sequential category change - can a transition happen at all
# Human Population and species introduction events
eventrule = let events=landscape_events.mus,
                states=states,
                history=history
                # D=dims(masks.mus),
    SetGrid{:landcover,:landcover}() do data, l1, l2
        current_year = currenttime(data)
        if hasselection(history, Ti(At(current_year)))
            foreach(eachindex(l1), l1, view(history, Ti(At(current_year)))) do I, state, hist
                # Fill water
                if hist.water && count(hist) == 1
                    l1[I] = states.water 
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
pressure = let preds=lc_predictions, ngridcells=size(sum(masks.mus))
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
    states,
    inertia=P(1.0),
    transitions,
    logic=(; direct=logic, indirect=indirect_logic(logic)),
    pressure,
    suitability=map(_ -> 1, states), #Aux{:suitability}(),
    history=Aux{:history}(),
    fixed=false,
    perturbation=P(2.0; bounds=(0.0, 10.0), label="perturbation"),
)
init_state = (; 
    landcover=Rasters.mask!(fill(1, dims(masks.mus); missingval=0), with=masks.mus)
)

aux = map(fix_order, (; history, suitability=suitability.mus))
tspan = 1600:2018
array_output = ArrayOutput(init_state;
    aux,
    tspan,
    store=false,
    mask=masks.mus,
    boundary=Remove(),
    padval=0,
)
ruleset = Ruleset(staterule, eventrule; proc=CPUGPU());
simdata = DynamicGrids.SimData(array_output, ruleset);

output = MakieOutput(init_state;
    aux, tspan,
    fps=100,
    store=true,
    mask=masks.mus,
    boundary=Remove(),
    padval=0,
    ruleset,
    sim_kw=(; printframe=true),
) do fig, frame, time
    axis1 = Axis(fig[1, 1])
    axis2 = Axis(fig[1, 2])
    linkaxes!(axis1, axis2)
    landcover = Observable(frame[].landcover)
    known_slices = Observable(view(striped_compiled, Ti(1)))
    on(frame) do f
        landcover[] = f.landcover
        t = tspan[time[]]::Int
        if hasselection(striped_compiled, Ti(At(t)))
            known_slices[] = view(striped_compiled, Ti(At(t)))
            notify(known_slices)
        end
        notify(landcover)
    end
    colormap = cgrad(:batlow, length(states)+1; categorical=true)
    hm = Makie.image!(axis1, landcover; colorrange=(first(states) -1.5, last(states) + 0.5), colormap, interpolate=false)
    Makie.image!(axis2, known_slices; colorrange=(first(states) -1.5, last(states) + 0.5), colormap, interpolate=false)
    ticks = (collect(0:length(states)), vcat(["mask"], collect(string.(propertynames(states)))))
    Colorbar(fig[1, 3], hm; ticks)
    return nothing
end

sim!(array_output, ruleset; simdata, proc=CPUGPU(), printframe=true);

predicted_lc = Rasters.combine(RasterSeries(array_output, dims(array_output)))
lc_predictions = map(NamedTuple(states)) do state
    rebuild(predicted_lc.landcover .== state; missingval=false, refdims=())
end |> RasterStack

lc_predictions_path = "$outputdir/lc_predictions.nc"
# write(lc_predictions_path, Rasters.modify(A -> UInt8.(A), lc_predictions))
lc_predictions = rebuild(Rasters.modify(BitArray, RasterStack(lc_predictions_path)); missingval=false)
# netcdf has the annoying center locus for time
lc_predictions = Rasters.set(lc_predictions, Ti => Int.(maybeshiftlocus(Start(), dims(lc_predictions, Ti), )))
