using CategoricalArrays
using GLM

categorical_soil = rebuild(soiltypes; data=categorical(parent(soiltypes)), name=:categorical_soiltypes)
soiltypes = norder_stack[:soiltypes]

# Flat island
# fi_selectors = X(57.64..57.69), Y(-19.9..(-19.86))
# plot(dems.mus[fi_selectors...])
# fi_pixels = length(collect(skipmissing(dems.mus[fi_selectors...])))
# step(parent(dims(dems.mus, X))) * 111
# fi_area = uconvert(u"km^2", npixels * 90.0u"m" * 90.0u"m") 
#
# Count pixels of land cover categories
cover_pixels = map(lc_masks) do island
    sort([map((k, v) -> k => sum(skipmissing(v)), keys(island), values(island))...]; by=last, rev=true)
end
cover_pixels.mus
cover_pixels.reu

# Fit land cover to soil, etc
using GLM
all_layers = dems, aspectrasters, sloperasters, distance_to_water

predictor_stacks = map(soilmasks, highway_masks, lc_masks, distance_to_ports, all_layers...) do hw, lc, dp, a...
    RasterStack(lc..., dp..., a...)
end;

st = RasterStack(norder_rasters..., veg_classes.mus..., predictor_stacks.mus..., categorical_soil)
df = DataFrame(st)

# model = lm(@formula(open_dry_palm_rich_woodland ~ rainfall + categorical_soiltypes + aspect + slope), df)
# r2(model)
# model = lm(@formula(semi_dry_evergreen_forest ~ DEM + categorical_soiltypes + ), df)
ag = rebuild(lc_masks.mus[:Sugarcane] .| lc_masks.mus[:Other_cropland]; name=:crops)
plot(ag; c=:viridis)
st = RasterStack(crops, soilmasks[:Latosolic_Red_Prairie_Soils], lc_masks.mus[:Barren_land], categorical_soil)

model = lm(@formula(Sugarcane ~ Latosolic_Red_Prairie_Soils), DataFrame(st))

model = lm(@formula(Barren_land ~ categorical_soiltypes), DataFrame(st))
model = lm(@formula(crops ~ categorical_soiltypes), DataFrame(st))
# model = lm(@formula(Discontinuous_urban ~ to_coasts), df)
# model = lm(@formula(Discontinuous_urban ~ to_coasts), df)
r2(model)
plot(soiltypes)

plot(soilmasks[:Latosolic_Red_Prairie_Soils]; c=:viridis)
plot(lc_masks.mus[:Discontinuous_urban]; c=:viridis)
plot(lc_masks.mus[:Discontinuous_urban] .& soilmasks[:Latosolic_Red_Prairie_Soils]; c=:viridis)
plot(crops .& soilmasks[:Latosolic_Red_Prairie_Soils]; c=:viridis)
fraction_used_urban = map(soilmasks) do A
    sum(skipmissing(lc_masks.mus[:Discontinuous_urban] .& A)) / sum(skipmissing(A))
end |> pairs |> collect
sort(fraction_used_urban; by=last, rev=true)
fraction_used_forest = map(soilmasks) do A
    sum(skipmissing((lc_masks.mus[:Forest]) .& A)) / sum(skipmissing(A))
end |> pairs |> collect
sort(fraction_used_forest; by=last, rev=true)
plot(lc_masks.mus[:Herbaceaous_vegetation] .| lc_masks.mus[:Discontinuous_urban]; c=:viridis)

fractions_by_lc = map(namedkeys(NamedTuple(lc_masks.mus))) do lckey
    s = sum(skipmissing(lc_masks.mus[lckey]))
    if s == 0
        nothing
    else
        fraction_used_ag = map(soilmasks) do A
            sum(skipmissing(lc_masks.mus[lckey] .& A)) / s
        end |> pairs |> collect |> x -> sort(x; by=x->x[2][1], rev=true)
    end
end |> pairs |> xs -> filter(x -> !isnothing(x[2]), xs)

fraction_by_soil = map(namedkeys(NamedTuple(soilmasks))) do soilkey
    map(lc_masks.mus) do A
        sum(skipmissing(soilmasks[soilkey] .& A)) / sum(skipmissing(soilmasks[soilkey]))
    end |> pairs |> collect |> x -> sort(x; by=x->x[2], rev=true)
end |> pairs

soilpixels = map(soilmasks) do A
    sum(skipmissing(A))
end |> pairs |> collect |> x -> sort(x; by=x->x[2], rev=true)

lcpixels = map(lc_masks.mus) do A
    sum(skipmissing(A))
end

sort_nt(x) = x |> pairs |> collect |> x -> sort(x; by=last, rev=true)


lc_masks.mus

lc_masks.mus[:Herbaceaous_vegetation] .| lc_masks.mus[:Discontinuous_urban]

plot(lc_masks.mus[:Herbaceaous_vegetation]; c=:viridis)
plot(lc_masks.mus[:Shrub_vegetation]; c=:viridis)
plot(lc_masks.mus[:Mangrove]; c=:viridis)

plot(lc_masks.mus[:Forest]; c=:viridis)
plot(lc_masks.mus[:Discontinuous_urban]; c=:viridis)
plot(lc_masks.mus[:Sugarcane]; c=:viridis)
plot(lc_masks.mus[:Other_cropland]; c=:viridis)
plot(lc_masks.mus[:Barren_land]; c=:viridis)

cats = lc_categories[(:Discontinuous_urban, :Forest, :Sugarcane, :Other_cropland, :Barren_land)]
plot(lc_rasterized.mus; c=cgrad(:Paired_12, 12; categorical=true), clims=(1, 12))


fractions_by_lc = map(namedkeys(NamedTuple(lc_masks2))) do lckey
    s = sum(skipmissing(lc_masks2[lckey]))
    if s < 2000
        nothing
    else
        fraction_used_ag = map(soilmasks) do A
            sum(skipmissing(lc_masks2[lckey] .& A)) / s
        end |> pairs |> collect |> x -> sort(x; by=x->x[2][1], rev=true)
    end
end |> pairs |> xs -> filter(x -> !isnothing(x[2]), xs)

fraction_by_soil = map(namedkeys(NamedTuple(soilmasks))) do soilkey
    x = map(lc_masks2) do A
        sum(skipmissing(soilmasks[soilkey] .& A)) / sum(skipmissing(soilmasks[soilkey]))
    end |> pairs |> collect |> x -> sort(x; by=x->x[2], rev=true)
end |> pairs

soiltypes1 = mask(norder_rasters[:soiltypes]; with=lc_rasterized.mus)
soilmasks1 = map(soilnums) do v
    soiltypes1 .== v
end |> RasterStack

lc_rasterized1 = mask(lc_rasterized.mus; with=soiltypes1)
lc_masks1 = map(namedkeys(lc_categories), lc_categories) do k, v
    Raster(mask(lc_rasterized1 .== v; with=soiltypes1, missingval=missing); name=k)
end
lc_masks1 = mask(mask(lc_masks1; with=soiltypes1); with=lc_rasterized.mus)

lc2 = clean_categories(lc_rasterized1;
    categories=cats, neighborhood=Moore{1}(), missingval=missingval(lc_rasterized1),
    keep_neigborless=true, despecle=false
)
lc2 = mask(lc2; with=soiltypes1)

# |> A -> classify(A, 
#    lc_categories[:Discontinuous_urban] => 1,
#    lc_categories[:Sugarcane] => 2,
#    lc_categories[:Other_cropland] => 2,
#    lc_categories[:Barren_land] => 0,
#    lc_categories[:Forest] => 1,
#    others=0
# )
plot(lc2; c=cgrad(:Paired_12, 12; categorical=true), clims=(1, 12))

lcpixels = map(RasterStack(lc_masks2...)) do A
    sum(skipmissing(A))
end |> sort_nt
lc_masks2 = map(namedkeys(lc_categories), lc_categories) do k, v
    Raster(mask(lc2 .== v; with=lc2, missingval=missing); name=k)
end

count(x -> true, skipmissing(soiltypes1))
count(x -> true, skipmissing(lc2))

o = [sum(skipmissing(lc .& s)) for lc in lc_masks2, s in soilmasks1]
e = [(sum(skipmissing(lc)) / sum(o)) * sum(skipmissing(s)) for lc in lc_masks2, s in soilmasks]
A = DimArray((o .- e) ./e, (lc=collect(keys(lc_masks2)), soil=collect(keys(soilmasks))))
Plots.heatmap(A, xrotation=30, clims=(-3, 3), c=:roma, size=(1000, 1000))

soil_lc_weight = [m[i, j] / sum(m) for i in axes(m, 1), j in axes(m, 2)]
sum(soil_lc)

watermasks.mus .& lc_masks.mus[:Sugarcane] |> plot


# cleared_rasters = map(lost_land_rasters) do i
#     cleared = maximum(i.rem)
#     (x -> x == cleared).(replace_missing(i.rem))
# end

# cleared_rem_plots = map(veg_phases, cleared_rasters) do vp, rem
#     plot(
#         plot((!).(vp[end]); c=:spring, title="Not cleared"), 
#         plot((!).(rem); c=:spring, title="Remnant")
#     )
# end

includet("cost_distance.jl")

m = typemin(Int)
costs = [
    1 3 4 4 3 2
    4 6 2 3 7 6         
    5 8 7 5 6 6
    1 4 5 m 5 1 
    4 7 5 m 2 6
    1 2 2 1 3 4
]
origins = [
    0 1 1 0 0 0
    0 0 1 0 0 0         
    0 0 0 0 0 0
    0 0 0 0 0 0 
    0 0 0 0 0 0
    2 0 0 0 0 0
]
cost_distance(; origins, costs, missingval=typemin(Int))

using ProfileView
using Cthulhu
using BenchmarkTools
@btime cost_distance(; origins, costs, missingval=typemin(Int))
@profview for _ in 1:10000 cost_distance(source, costs; missingval=typemin(Int)) end

