using LandscapeChange
using Distributions
using DimensionalData
using DimensionalData.LookupArrays
using StatsPlots
using Unitful
using StaticArrays
using ReverseStackTraces
using GLMakie
using Statistics, StatsBase, GLM
using ThreadsX

# From woods and forests in Mauritius, p 40:40
# 1880: 70,000 acres out of 300,000 remain
# 35,000 were native (but "dilatipated and ruined?")
# Also mentions that invasives replace natives
includet("raster_common.jl")
includet("roads.jl")
include("travel_cost.jl")
include("tabular_data.jl")
include("/home/raf/.julia/dev/LandscapeChange/src/makie.jl")

states = NamedVector(lc_categories)

# homiisland = Raster("/home/raf/PhD/Mascarenes/Data/Generated/Landcover/mus_landcover.tif")
# homiisland2 = resample(homiisland; to=dems.mus)
# Rasters.rplot(homiisland2 .== 8)
# write("/home/raf/PhD/Mascarenes/Data/Generated/Landcover/mus_landcover_2.tif")

println("Getting timeline slices...")
# includet("map_file_list.jl")
includet("map_file_list_2.jl")
files = define_map_files()
slices = make_raster_slices(masks, lc_categories; category_names);
# Rasters.rplot(slices.mus.timelines.abandoned; colorrange=(0, 1))
nv_ser = namedvector_raster.(slices.mus.timeline)
# sandwich!(nv_ser)
striped = stripe_raster(nv_ser, states)
# Rasters.rplot(striped[5:end]; colorrange=(1, 6))

const P = Param 
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
# to category from category with intermediate steps
indirect_logic = map(states) do s1
    map(states) do s2
        can_change(states, logic, s1, s2)
    end
end
pairs(indirect_logic)

slices = make_raster_slices(masks, lc_categories; category_names);
# Rasters.rplot(slices.mus.timelines.abandoned; colorrange=(0, 1))
nv_ser = namedvector_raster.(slices.mus.timeline)
striped = stripe_raster(nv_ser, states)
compiled = compile_timeline(indirect_logic, nv_ser; fill=true, continuity=false, transition=false)
compiled = compile_timeline(indirect_logic, nv_ser; fill=true, continuity=false, transition=false)
simplified = compile_timeline(indirect_logic, nv_ser; fill=true, continuity=true, transition=false)
trans = compile_timeline(indirect_logic, nv_ser; fill=true, continuity=false, transition=true)
complete = compile_timeline(logic, nv_ser; fill=true, continuity=true, transition=true)
striped_compiled = stripe_raster(compiled.timeline, states)
striped_simplified = stripe_raster(simplified.timeline, states)
striped_complete = stripe_raster(complete.timeline, states)
striped_trans = stripe_raster(trans.timeline, states)
striped_error = stripe_raster(complete.error, states)
# diff = map((x, y) -> (x .& .!(y)), simplified.timeline, compiled.timeline)
# striped_diff = stripe_raster(diff, states)
rnge = Ti(4:size(striped, Ti))
colorrange = (1, 6)
fig = Figure()
empty!(fig); Rasters.rplot(fig[1, 1], striped[rnge]; colorrange)
empty!(fig); Rasters.rplot(fig[1, 1], striped_compiled[rnge]; colorrange)
empty!(fig); Rasters.rplot(fig[1, 1], striped_simplified[rnge]; colorrange)
empty!(fig); Rasters.rplot(fig[1, 1], striped_trans[rnge]; colorrange)
empty!(fig); Rasters.rplot(fig[1, 1], striped_complete[rnge]; colorrange)
empty!(fig); Rasters.rplot(fig[1, 1], striped_error[rnge]; colorrange)

a = NV(native=false, cleared=false, abandoned=false, urban=true, forestry=true, water=true)
b = NV(native=false, cleared=true, abandoned=false, urban=false, forestry=false, water=false)
_merge_all_possible(a, b, logic) |> pairs
_merge_all_possible(b, a, logic) |> pairs

tbl = Shapefile.Table("/home/raf/PhD/Mascarenes/Data/Claudia/Black River Gorges and other shapes/npcs.shp") |> DataFrame
tbl = Shapefile.Table("/home/raf/PhD/Mascarenes/Data/Claudia/Demo/GIS WILD LIFE FOUNDATION/SHAPE FILE/land_use_WGS_region.shp") |> DataFrame
union(tbl.LAND_USE)
swamps = filter(:LAND_USE => x -> !ismissing(x) && x == "Marsh or Swamp", tbl)
cleared = filter( :LAND_USE => 
  x -> !ismissing(x) && x in ("Other Plantation", "Sugar Cane", "Other PlanTation", "Tea Plantation", "Sugar cane"), tbl)
# swamps_rast = boolmask(swamps; to=masks.mus)
cleared_rast = boolmask(cleared; to=masks.mus)
other_rast = .!(cleared_rast)
lc = cleared_rast .* 1 .+ other_rast .* 2 .* masks.mus
write("/home/raf/PhD/Mascarenes/Data/Generated/Landcover/mus_wlf_shape.tif", lc)



# for pdf in readdir("/home/raf/PhD/Mascarenes/maps/Mauritius/Studies_of_Mascarine_birds"; join=true)
#     name, ext = splitext(pdf)
#     ext == ".pdf" || continue
#     png = name * ".png" 
#     run(`pdftoppm $pdf $png -png -r 300`)
# end

function countcats(data, categories=union(data))
    map(categories) do cat
        cat => count(==(cat), data)
    end
end
Plots.plot!(human_pop_timeline.reu)
Plots.plot(human_pop_timeline.mus)
human_pop_timeline.mus

# best_slices = lc[At([1600, 1709, 1723, 1772, 1810, 1905, 1992])]
# Plots.plot(best_slices)
DimensionalData.bounds(lookup(slices.mus.timeline, Ti))
eltype(slices.mus.timeline)
cats1 = (:cleared, :abandoned, :urban, :forestry, :water)
cat_counts = map(slice(complete.timeline, Ti)) do slice
    total_counts = zeros(Int, size(first(slice)))
    known_counts = zeros(Int, size(first(slice)))
    for categories in slice
        total_counts .+= categories
        if count(categories) == 1
            known_counts .+= categories
        end
    end
    return total_counts, known_counts
end

map(first.(cat_counts), last.(cat_counts)) do t, k
    collect(zip(keys(states), k, t))
end

Plots.scatter(human_pop_timeline.mus)
pop_cat_counts = map(slices.mus.timelines[cats1], cat_counts) do ct, cc
    map(lookup(ct, Ti)) do t
        (; count=cc[Contains(t)], pop=human_pop_timeline.mus[At(t)])
    end
end

# Cleared land is used for urbanisation by 1992, so don't use it in the model
cleared_model = lm(@formula(count ~ pop^2 + pop), pop_cat_counts.cleared)
urban_model = lm(@formula(count ~ pop^2), pop_cat_counts.urban)
ti = dims(human_pop_timeline.mus, Ti)
pops = map(pop -> (; pop), human_pop_timeline.mus)
Plots.plot(first.(pops))

cleared_pred = DimArray(predict(cleared_model, parent(pops)), ti)
urban_pred = DimArray(predict(urban_model, parent(pops)), ti)
lc_predictions = map((cleared, urban) -> (; cleared, urban), cleared_pred, urban_pred)
Plots.plot(cleared_pred)
Plots.scatter!(cat_counts.urban)
Plots.plot!(urban_pred)
Plots.scatter!(cat_counts.cleared)
Plots.plot!(human_pop_timeline.mus)

b = (; bounds=(0.0, 2.0))
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
fixed = false # or a raster mask of fixed landcover
human_suitability = map(travel_times, slope_stacks, dems) do tt, ss, dem
    slope_suitability = (1 .- ss.slope) .^ 2 # Steep bad, flat good
    slices = map(tt) do travel_time
        travel_suitability = (1 .- travel_time / 30u"hr")
        reverse(replace_missing(travel_suitability .* slope_suitability, 0.0); dims=Y)
    end
    Rasters.combine(slices, Ti)
end
# Plots.plot(travel_times.rod)

distance_to_water = map(island_keys) do k
    fix_order(Raster(joinpath(distancedir, string(k), "to_water.tif")))
end
# plot(suitability.mus; clims=(0, 1), legend=false)
# plot(suitability.reu; clims=(0, 1), legend=false)
suitability = map(human_suitability, distance_to_water) do hs, dtw
    # dtw = 1 ./ (1 .+ sqrt.(replace_missing(dtw, Inf)))
    native = fill!(similar(hs), 1.0)
    cleared = hs .^ 2# .* dtw
    abandoned = 1 .- hs
    urban = hs .^ 21# .* dtw
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
#
suit = human_suitability.rod
pop = pop_density.rod
pop_models = map(human_suitability, pop_density) do suit, pop
    ag_pop_density = Rasters.aggregate(sum, pop, 20; skipmissingval=true)
    ag_pop_density = replace(ag_pop_density, NaN=>0.0)
    resample_suit = mask!(resample(suit[Ti=1]; to=ag_pop_density); with=ag_pop_density)
    mask!(ag_pop_density; with=resample_suit)
    df = RasterStack((pop=replace_missing(ag_pop_density), suit=replace_missing(resample_suit)))
    model = lm(@formula(pop ~ suit), df)
end
map(r2, pop_models)
# Plots.plot(dems.mus)

# Rasters.aggregate(sum, pop_density.reu, 50; skipmissingval=true) |> Plots.plot



# The amount of influence neighbors have
b = (; bounds=(1.00, 5.0))
transitions = (
    native = (
        native=0.0,#Exponential(P(1.7; b..., label="native from native")),
        cleared=0.0,#Exponential(P(1.0; b...)),
        abandoned=0.0,
        urban=0.0,
        forestry=0.0,#Exponential(P(1.0; b...)),
        water=0.0,
    ),
    cleared = (
        native=0.0,
        cleared=Exponential(P(1.7; b..., label="cleared from cleared")),
        abandoned=0.0,#Exponential(P(1.0; b...)),
        urban=Exponential(P(1.0; b..., label="cleared from urban")),
        forestry=0.0,#Exponential(P(1.0; b...)),
        water=0.0,
    ),
    abandoned = (
        native=0.0,
        cleared=0.0,#Exponential(P(1.0; b...)),
        abandoned=Exponential(P(1.7; b..., label="abandoned from abandoned")),
        urban=0.0,#Exponential(P(1.0; b...)),
        forestry=0.0,#Exponential(P(1.0; b...)),
        water=0.0,
    ),
    urban = (
        native=0.0,
        cleared=Exponential(P(1.7; b..., label="urban from cleared")),
        abandoned=0.0,#Exponential(P(1.0; b...)),
        urban=Exponential(P(1.7; b..., label="urban from urban")),
        forestry=0.0,#Exponential(P(1.0; b...)),
        water=0.0,
    ),
    forestry = (
        native=0.0,
        cleared=0.0,
        abandoned=0.0,
        urban=0.0,
        forestry=Exponential(P(1.7; b..., label="forestry from forestry")),
        water=0.0,
    ),
    water = (
        native=0.0,
        cleared=0.0,
        abandoned=0.0,
        urban=0.0,
        forestry=0.0,
        water=0.0,
    ),
)

# The logic of sequential category change - can a transition happen at all
# Human Population and species introduction events
# eventrule = let events=landscape_events.mus,
#                 states=states,
#                 D=dims(masks.mus),
#                 history=history,
#                 lc=lc,
#     SetGrid{:landcover,:landcover}() do data, l1, l2
#         current_year = currenttime(data)
#         foreach(history, ) do timeseries
#             l2 .= lc[At(current_year)]
#         end
#         # for event in events
#         #     if event.year == current_year
#         #         p = event.geometry
#         #         I = DimensionalData.dims2indices(D, (X(Contains(p.X)), Y(Contains(p.Y))))
#         #         Iu = map(i -> i-1:i+1, I)
#         #         Ic = map(i -> i-10:i+10, I)
#         #         l1[Iu...] .= states.urban
#         #         l1[Ic...] .= states.cleared
#         #         l2[Iu...] .= states.urban
#         #         l2[Ic...] .= states.cleared
#         #         # l2 .= states.cleared
#         #     end
#         # end
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
    leverage=P(1.0; bounds=(0.0, 20.0))
    smoothing=P(1.0; bounds=(1.0, 5.0))
    # cleared=P(1.5; b...),
    abandoned=P(1.0; b...)
    # urban=P(1.4; b...)
    forestry=P(1.0; b...)
    native = 0.0
    water = 0.0
    (data, rule) -> begin
        predicted = preds[At(currenttime(data))]
        stategrid = data[:landcover]
        # ngridcells = if isnothing(DynamicGrids.mask(data))
        #     length(stategrid)
        # else
        #     sum(DynamicGrids.mask(data))
        # end
        ncleared = ThreadsX.sum(==(rule.states.cleared), parent(stategrid); init=0.0)
        nurban = ThreadsX.sum(==(rule.states.urban), parent(stategrid); init=0.0)
        # println(stdout, (; ncleared, pcleared=predicted.cleared, nurban, purban=predicted.urban))
        urban = leverage * sign(predicted.urban - nurban) * 
            max(0.0, predicted.urban) / max(nurban, max(0.0, predicted.urban) / smoothing)
        cleared = leverage * sign(predicted.cleared - ncleared) * 
            max(0.0, predicted.cleared) / max(ncleared, max(0.0, predicted.cleared) / smoothing)
        NamedVector(; native, cleared, abandoned, urban, forestry, water)
    end
end

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
    logic=(; direct=logic, indirect=indirect_logic),
    pressure,
    suitability=Aux{:suitability}(),
    fixed=false,
    perturbation=P(0.1; bounds=(0.0, 1.0), label="perturbation"),
)
staterule
init_state = (; 
    landcover=Rasters.mask!(fill(1, dims(masks.mus); missingval=0), with=masks.mus)
)
cats1 = (:cleared, :abandoned, :urban, :forestry, :water)
cats1_nt = NamedTuple{cats1}(cats1)
native_timedim = Ti(rebuild(lookup(slices.mus.timelines[:cleared], Ti); data=[1600, 2000]))
mus_native_1600 = masks.mus
mus_native_1999 = set(fix_order(boolmask(mus_native_veg_1999)), X=>dims(mus_native_1600, X), Y=>dims(mus_native_1600, Y))
native_history = cat(mus_native_1600, mus_native_1999; dims=native_timedim)
history = merge((; native=native_history), map(cats1_nt) do k
    Rasters.combine(slices.mus.timelines[k])
end)

map(missingval, history)
map(eltype, history)
# Rasters.rplot(history; colorrange=(1, 5))
aux = map(fix_order, (; suitability=suitability.mus, history...))
result_output = ResultOutput(init_state;
    aux,
    tspan=1630:1:1931,
    store=false,
    mask=masks.mus,
    boundary=Remove(),
    padval=0,
)
ruleset = Ruleset(staterule; proc=CPUGPU());
simdata = DynamicGrids.SimData(result_output, ruleset);
@time sim!(result_output, ruleset; simdata, proc=CPUGPU()); 

# using ProfileView
# using Cthulhu
# @profview for i in 1:100000 DynamicGrids.descendable(simdata) end
# @descend DynamicGrids.descendable(simdata)
# @profview sim!(result_output, ruleset; simdata, proc=CPUGPU()); 

# plot(RasterStack(slices.mus.timelines.urban); legend=false, size=(1000, 850), margin=-1mm)
# plot(RasterStack(slices.mus.timelines.lc); legend=false, clims=(0, 5), size=(1000, 850), margin=-1mm)
# plot(last(slices.mus.timelines.lc); legend=false, clims=(0, 5), size=(1000, 850), margin=-1mm)
# savefig("mus_landcover_maps.png")

# lc_slices = collect(slices.mus.timelines.lc)
# plot(last(lc_slices))
# lc_slices
# uncertain_lc = mask(last(lc_slices); with=replace_missing(lc_2017.mus.uncertain, false))
# plot(uncertain_lc)
# countcats(uncertain_lc, lc_categories) |> pairs
# slices.mus.timelines
# lc_years = map(keys(slices.mus.timelines.lc)) do key
#     parse(Int, string(key)[end-3:end])
# end |> collect
# lc_timeline = RasterSeries(lc_slices, Ti(lc_years))
# lc_fractions = map(A -> LandscapeChange.cover_fraction(A; categories=lc_categories), lc_timeline)

# Plots.plot(getproperty.(lc_fractions, :abandoned); labels="abandoned")
# Plots.plot!(getproperty.(lc_fractions, :cleared); labels="cleared")
# Plots.plot!(getproperty.(lc_fractions, :native); labels="native")
# Plots.plot!(getproperty.(lc_fractions, :forestry); labels="forestry")
# Plots.plot!(getproperty.(lc_fractions, :urban); legend=:left, labels="urban", title="Mauritius land cover history")
# Plots.plot!(twinx(), human_pop_timeline.mus; legend=:top, color=:black, labels="human population")
# savefig("mus_landcover_hist.png")

# pop_at_slice = human_pop_timeline.mus[At(lookup(lc_fractions, Ti))]

# table = map(lc_fractions, pop_at_slice) do f, pop
#     NamedTuple(merge(f, (; pop)))
# end |> DataFrame

# f = @formula(x ~ a + b + c + d)
# length(pop_at_slice)
# model = lm(@formula(1 - native ~ pop^2 + pop), table)
# @formula(1 - native ~ pop^2 + pop)
# pops = DataFrame((pop = 1:100:1000000,))
# plot(pops.pop, predict(model, pops))
# Plots.scatter!(table.pop, 1 .- table.native)

# persistence = map(GeometryBasics.TupleView{2,1}(lc_slices)) do xs
#     LandscapeChange.cover_persistence(xs...; categories=lc_categories)
# end
# persistence_timeline =
#     DimArray(persistence, set(dims(lc_timeline, Ti)[begin:end-1], Sampled(; sampling=Intervals(Start()), span=Irregular(first(lc_years), last(lc_years)))))
# annual_persistence = map(enumerate(persistence_timeline)) do (i, t)
#     b = cellbounds(persistence_timeline, i)[1]
#     years = b[2] - b[1]
#     1 .- ((1 .- t) ./ years)
# end
# annual_persistence_timeline = rebuild(persistence_timeline, annual_persistence)

# Plots.plot(getproperty.(annual_persistence_timeline, :abandoned); labels="abandoned")
# Plots.plot!(getproperty.(annual_persistence_timeline, :cleared); labels="cleared")
# Plots.plot!(getproperty.(annual_persistence_timeline, :native); labels="native")
# Plots.plot!(getproperty.(annual_persistence_timeline, :urban); labels="urban", title="Annual persistence of land cover classes", ylims=(0.96, 1.0))
# Plots.plot!(getproperty.(annual_persistence_timeline, :forestry); legend=:left, labels="forestry", title="Annual persistence of land cover classes", ylims=(0.96, 1.0))
# Plots.plot!(twinx(), human_pop_timeline.mus; legend=:bottomleft, color=:black, labels="human population")
# savefig("annual_persistence.png")

# transition_probs = map(GeometryBasics.TupleView{2,1}(lc_slices)) do xs
#     LandscapeChange.cover_change(xs...; categories=lc_categories)
# end
# transition_timeline =
#     DimArray(transition_probs, set(dims(lc_timeline, Ti)[begin:end-1], Sampled(; sampling=Intervals(Start()), span=Irregular(first(lc_years), last(lc_years)))))
# annual_transition = map(enumerate(transition_timeline)) do (i, ts)
#     b = cellbounds(transition_timeline, i)[1]
#     years = b[2] - b[1]
#     map(ts) do t
#         t ./ years
#     end
# end
# annual_transition_timeline = rebuild(transition_timeline, annual_transition)
# annual_transition_timeline[1]
# getproperty.(annual_transition_timeline, :native)
# x = annual_transition_timeline[1]


Rasters.rplot(history.abandoned)
Rasters.rplot(history.cleared)
Rasters.rplot(history.urban)
Rasters.rplot(history.native)
Rasters.rplot(history.forestry)

Rasters.rplot(slices.mus.files.landcover_1965.grouped.cleared)


# set_theme!(theme_dark())
#
# ps = map(lookup(lc, Ti)) do t
#     r = replace_missing(lc[Ti(At(t))], NaN)
#     figure = Figure(; 
#         backgroundcolor=:transparent,
#     )
#     axis = Axis(figure[1, 1];
#         xticklabelsvisible=false, 
#         yticklabelsvisible=false,
#         xticksvisible=false, 
#         yticksvisible=false,
#         xgridvisible=false,
#         ygridvisible=false,
#         bottomspinevisible=false,
#         topspinevisible=false,
#         leftspinevisible=false,
#         rightspinevisible=false,
#         aspect=DataAspect(),
#     )
#     text = "$t"
#     textpos = 57.3, -20.0
#     Makie.text!(axis, textpos...; text, fontsize=80, align=(:left, :top)) 
#     p = Makie.heatmap!(axis, r; 
#         colormap=:batlow, 
#         colorrange=(0, 6),
#     )jj\zz
#     save("images/mus_lc_$t.png", figure)
#     p
# end

# CairoMakie.activate!()
# GLMakie.activate!()
# savefig("mauritius_timeline.png")
# Rasters.rplot(Rasters.combine(slices.mus.timelines.lc, Ti); 
#     colormap=:batlow, colorrange=(0, 5),
#     aspect=DataAspect(),
#     axis = (
#         xticklabelsvisible=false, 
#         yticklabelsvisible=false,
#         xticksvisible=false, 
#         yticksvisible=false,
#         xgridvisible=false,
#         ygridvisible=false,
#         bottomspinevisible=false,
#         topspinevisible=false,
#         leftspinevisible=false,
#         rightspinevisible=false,
#     ),
# )
# lc = map(slices.mus.timelines.lc) do A
#     reverse(A; dims=Y)
# end |> x -> set(x, Ti=>Intervals(End())) |> x -> set(x, Ti=>Irregular((0, 1992)))
# a = @animate for A in lc
#     Plots.plot(A; legend=false, clims=(1, 5))
# end
# Plots.gif(a, "timeseries.gif", fps=1)
# Plots.plot(lc; size=(2000,1300), legend=false, clims=(1, 5))
# savefig("mauritius_timeline.png")

