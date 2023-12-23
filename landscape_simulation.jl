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
using GLMakie
using Rasters.LookupArrays

# From woods and forests in Mauritius, p 40:40
# 1880: 70,000 acres out of 300,000 remain
# 35,000 were native (but "dilatipated and ruined?")
# Also mentions that invasives replace natives
include("travel_cost.jl")
include("tabular_data.jl")
include("map_file_functions.jl")

states = NamedVector(lc_categories)

const P = RealParam
const NV = NamedVector

# to category from category
logic = NV(
    native    = NV(native=true,  cleared=false, abandoned=false, urban=false, forestry=false,  water=false),
    cleared   = NV(native=true,  cleared=true,  abandoned=true,  urban=false, forestry=false,  water=false),
    abandoned = NV(native=false, cleared=true,  abandoned=true,  urban=false,  forestry=false,  water=false),
    urban     = NV(native=true,  cleared=true,  abandoned=true,  urban=true,  forestry=false,  water=false),
    forestry  = NV(native=true,  cleared=true,  abandoned=true,  urban=false, forestry=true,   water=false),
    water     = NV(native=true,  cleared=true,  abandoned=true,  urban=true,  forestry=false,  water=true),
)

include("map_file_list_2.jl")
slices = compile_timeline(define_map_files(), masks, keys(lc_categories))
force = NV{propertynames(logic)}(map(x -> x in (:cleared, :urban, :water), propertynames(logic)))
nv_rasts = map(slices) do island
    Rasters.combine(namedvector_raster.(island.timeline))
end
striped_raw = map(nv_rasts) do nv_rast
    stripe_raster(nv_rast, states)
end
compiled = map(nv_rasts) do nv_rast
    cross_validate_timeline(logic, nv_rast; simplify=true, cull=true)#, force)
end
striped_compiled = map(compiled) do island
    stripe_raster(island.timeline, states)
end
striped_error = map(compiled) do island
    stripe_raster(island.error, states)
end

# Rasters.rplot(striped_raw.mus; colorrange=(1, 6))
# Rasters.rplot(striped_compiled.mus[Ti=1:16]; colorrange=(1, 6))
# Rasters.rplot(striped_raw.reu; colorrange=(1, 6))
# Rasters.rplot(striped_compiled.reu; colorrange=(1, 6))
# Rasters.rplot(striped_raw.rod; colorrange=(1, 6))
# Rasters.rplot(striped_compiled.rod; colorrange=(1, 6))
# b = map(x -> x.cleared, compiled.rod.timeline)
# Rasters.rplot(b)

cat_counts = let states=states
    map(human_pop_timelines, compiled) do human_pop, history
        map(slice(history.timeline[Ti=1600..2017], Ti)) do slice
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
                pop = human_pop[Near(year)]
                (; total, known, ratio=known/total, meancount=(known + total) / 2, year, pop)
            end
            nv = NamedVector{propertynames(states)}(vals)
        end
    end
end

high_certainty = map(cat_counts) do cc
    map(category_names) do k
        vals = map(cc) do val
            val[k]
        end
        filter(v -> v.ratio > 0.3, vals)
    end
end

# Cleared land is used for urbanisation by 1992, so don't use it in the model
lc_targets = map(high_certainty[(:mus,)]) do hc
    cleared_model = lm(@formula(meancount ~ pop^2 + pop), DataFrame(hc.cleared))
    urban_model = lm(@formula(meancount ~ pop^2), DataFrame(hc.urban))
    ti = dims(human_pop_timelines.mus, Ti)
    pops = map(pop -> (; pop), human_pop_timelines.mus)
    cleared_pred = DimArray(predict(cleared_model, parent(pops)), ti)
    urban_pred = DimArray(predict(urban_model, parent(pops)), ti)
    map((cleared, urban) -> (; cleared, urban), cleared_pred, urban_pred)
end
# Plots.plot(cleared_pred)
# Plots.scatter!(getproperty.(high_certainty.cleared, :known))
# Plots.scatter!(getproperty.(high_certainty.cleared, :meancount))
# Plots.scatter!(getproperty.(high_certainty.cleared, :total))
# Plots.plot(urban_pred)
# Plots.scatter!(getproperty.(high_certainty.urban, :known))
# Plots.scatter!(getproperty.(high_certainty.urban, :meancount))
# Plots.scatter!(getproperty.(high_certainty.urban, :total))
# Plots.plot!(human_pop_timelines.mus)

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
human_suitabilities = map(travel_times, slope_stacks, dems) do tt, ss, dem
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

suitabilities = map(human_suitabilities, distance_to_water) do hs, dtw
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

# mus_native_veg_tif_path = "/home/raf/PhD/Mascarenes/Data/Generated/mus_native_veg.tif"
# target_native_fraction = Raster(mus_native_veg_tif_path) ./ 4

include("landscape_rules.jl");

init_states = map(masks) do mask
    (;  landcover=Rasters.mask!(fill(1, dims(mask); missingval=0), with=mask),
        native_fraction=rebuild(fill(1.0, dims(mask)); missingval=nothing)
    )
end
auxs = map(landscape_events, compiled, suitabilities) do events, comp, suitability
    (; events, map(fix_order, (; history=comp.timeline, suitability))...)#, target_native_fraction))
end
tspans = (mus=1600:2020, reu=1600:2020, rod=1600:2020)
output_kw = map(init_states, masks, auxs, tspans) do init, mask, aux, tspan
    (; aux, mask, tspan, store=false, boundary=Remove(), padval=0)
end
array_outputs = map(init_states, output_kw) do init, kw
    ArrayOutput(init; kw...)
end
foreach(array_outputs) do output
    sim!(output, ruleset; printframe=true);
end

lc_predictions = map(array_outputs) do array_output
    predicted_lc = Rasters.combine(RasterSeries(array_output, dims(array_output)))
    map(NamedTuple(states)) do state
        rebuild(predicted_lc.landcover .== state; missingval=false, refdims=())
    end |> RasterStack
end

# include("raster_common.jl")
# foreach(lc_predictions, island_keys) do lc, k
#     lc_stack_path = "$outputdir/lc_predictions_$k.nc"
#     write(lc_stack_path, Rasters.modify(A -> UInt8.(A), lc))
# end
lc_predictions = map(island_keys) do k
    lc_stack_path = "$outputdir/lc_predictions_$k.nc"
    st = rebuild(Rasters.modify(BitArray, RasterStack(lc_stack_path)); missingval=false)
    Rasters.set(st, Ti => Int.(maybeshiftlocus(Start(), dims(st, Ti), )))
end
# netcdf has the annoying center locus for time

k = :reu
output = MakieOutput(getproperty(init_states, k);
    getproperty(output_kw, k)...,
    fps=100, store=false,
    ruleset, sim_kw=(; printframe=true),
) do (; layout, frame, time)
    axis1 = Axis(layout[1, 1])
    axis2 = Axis(layout[1, 2])
    axis3 = Axis(layout[1, 3])
    linkaxes!(axis1, axis2)
    landcover = Observable(frame[].landcover)
    native_fraction = Observable(frame[].native_fraction)
    known_slices = Observable(view(getproperty(striped_compiled, k), Ti(1)))
    on(frame) do f
        landcover[] = f.landcover
        native_fraction[] = f.native_fraction
        t = tspans.rod[time[]]::Int
        if hasselection(striped_compiled, Ti(At(t)))
            known_slices[] = view(striped_compiled.mus, Ti(At(t)))
            notify(known_slices)
        end
        notify(landcover)
        notify(native_fraction)
    end
    colormap = cgrad(:batlow, length(states)+1; categorical=true)
    hm = Makie.image!(axis1, landcover; colorrange=(first(states) -1.5, last(states) + 0.5), colormap, interpolate=false)
    nf = Makie.image!(axis2, native_fraction; colorrange=(0, 1), interpolate=false)
    Makie.image!(axis3, known_slices; colorrange=(first(states) -1.5, last(states) + 0.5), colormap, interpolate=false)
    # ticks = (collect(0:length(states)), vcat(["mask"], collect(string.(propertynames(states)))))
    # Colorbar(fig[1, 3], hm; ticks)
    return nothing
end




mus_veg_path = "/home/raf/PhD/Mascarenes/Data/Selected/Mauritius/Undigitised/page33_mauritius_vegetation_colored.tif"
reu_veg_path = "/home/raf/PhD/Mascarenes/Data/Dominique/Vegetation_Rasters/pastveg3.tif"
original_veg = (;
    mus=replace_missing(Raster(mus_veg_path), 0),
    reu=resample(replace_missing(Raster(reu_veg_path), 0); to=masks.reu),
)

using DataFrames
nf_slices = getproperty.(output, :native_fraction);
Rasters.combine(RasterSeries(nf_slices, dims(nf_slices)), Ti);
veg_change = map(lc_predictions[keys(original_veg)], original_veg) do p, v
    rebuild(UInt8.(broadcast_dims(*, p.native, v)); missingval=0)
end
# mus_veg_change = rebuild(UInt8.(broadcast_dims(*, lc_predictions.mus.native, mus_veg)); missingval=0)

Rasters.rplot(veg_change)
habitat_names = ["semi-dry_evergreen_forest", "open_dry_palm-rich_woodland", "wet_forest", "pandanus_swamp", "mossy_rainforest", "mangrove", "wetland vegetation"]
length(habitat_names)
habitat_sums = cat(map(1:7) do habitat
    dropdims(sum(==(habitat), veg_change; dims=(X, Y)); dims=(X, Y))
end...; dims=Dim{:habitat}(habitat_names))
init_total = sum(habitat_sums) do x
    x[1]
end
cum = cumsum(habitat_sums; dims=2)
x = lookup(habitat_sums, Ti)
fig = Figure()
ax = Axis(fig[1, 1])
for i in 7:-1:1
    y = parent(cum[habitat=i])
    Makie.lines!(x, y; color=:black)
    band!(x, fill(0, length(x)), y; label = "Label")
end
fig[1, 2] = Legend(fig, ax, habitat_names)



# timeline = NV{(:native, :cleared, :abandoned, :urban, :forestry, :water)}.[
#     (true, false, false, false, false, false)
#     (false, true, false, false, false, false)
#     (false, false, false, true, false, false)
#     (false, true, true, false, false, false)
#     (false, false, false, true, false, false)
# ]
# removed = NV{keys(states)}.([
#     (true, false, false, false, false, false)
#     (false, true, false, false, false, false)
#     (false, true, false, true, false, false)
#     (false, true, false, true, false, false)
#     (false, false, false, true, false, false)
# ])
# corrected = NV{keys(states)}.([
#     (true, false, false, false, false, false)
#     (false, true, false, false, false, false)
#     (false, true, false, true, false, false)
#     (false, false, false, true, false, false)
#     (false, false, false, true, false, false)
# ])

# transitions = logic
# reversed = LandscapeChange.reverse_transitions(transitions)
# indirect = LandscapeChange.indirect_transitions(transitions)
# reversed_indirect = LandscapeChange.reverse_transitions(indirect)
# force = NV{k}(map(x -> x in (:native, :urban, :water), propertynames(transitions)))
# rev = reverse(eachindex(timeline))

# timeline = vec(nv_rasts.rod[X=Near(63.393), Y=Near(-19.7392)])
# timeline1 = copy(timeline); LandscapeChange._apply_transitions!(timeline1, reversed, reversed_indirect, force)

# timeline2 = copy(timeline); LandscapeChange._apply_transitions!(view(timeline2, rev), transitions, indirect, force)
# timeline
# timeline1
# timeline2
# LandscapeChange._remove_intermediate_uncertainty!(timeline1)
# LandscapeChange._remove_intermediate_uncertainty!(timeline2)
# possible = map(.|, timeline1, timeline2)
# t, r = LandscapeChange.cross_validate_timeline!(copy(timeline), transitions; simplify=true, cull=false); t
# t, r = LandscapeChange.cross_validate_timeline!(copy(timeline), transitions; simplify=true, cull=true); t

# a = NV{propertynames(transitions)}((false, false, true, false, true, false))
# b = NV{propertynames(transitions)}((true, false, true, false, true, false))

# LandscapeChange._merge_all(a, b, reversed, reversed_indirect, force)
