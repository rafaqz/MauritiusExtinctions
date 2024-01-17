using ModelParameters
using StaticArrays
using DynamicGrids
using Dispersal
using LandscapeChange
using CSV
using DataFrames
using XLSX
using TerminalPager
using Rasters
using GLMakie
using NCDatasets
using Geomorphometry
using Stencils
using Statistics

# includet("optimisation.jl")
include("species_rules.jl")
include("raster_common.jl")
include("makie.jl")

function gpu_cleanup(A)
    x, y, ti = lookup(A, (X, Y, Ti))
    Rasters.set(A, 
        Ti => Sampled(first(ti) - 0.5:last(ti) - 0.5; sampling=Intervals(Start())),
        X => LinRange(first(x), last(x), length(x)), 
        Y => LinRange(first(y), last(y), length(y)),
    )
end
uncleared = gpu_cleanup(Rasters.modify(BitArray, Raster("uncleared.nc")))
# forested = gpu_cleanup(modify(BitArray, Raster("forested.nc")))

pred_df = CSV.read("tables/animals.csv", DataFrame)
introductions_df = CSV.read("tables/introductions.csv", DataFrame)
mascarene_species_csv = "tables/mascarene_species.csv"
# @async run(`libreoffice $mascarene_species_csv`)
all_species = CSV.read(mascarene_species_csv, DataFrame) |> 
    x -> subset(x, :Species => ByRow(!ismissing))
island_tables = map((mus=:mus, reu=:reu, rod=:rod)) do key
    df = DataFrame(subset(all_species, key => x -> .!ismissing.(x)))
    df.presence = df[!, key]
    df.extinct = map(df[!, "$(key)_extinct"]) do e
        ismissing(e) ? missing : eval(Meta.parse(e))::UnitRange
    end
    df.introduced = map(df[!, "$(key)_introduced"]) do e
        ismissing(e) ? missing : eval(Meta.parse(e))::UnitRange
    end
    df
end

island_extinct_tables = map(island_tables) do tbl
    DataFrame(subset(tbl, :extinct => x -> .!ismissing.(x)))
end

get_species_names(table) = Tuple(Symbol.(replace.(skipmissing(table.Species), Ref(" " => "_"))))
island_names = NamedTuple{keys(island_tables)}(keys(island_tables))
island_extinct = map(island_tables, island_names) do table, name
    subset(table, Symbol(name, "_extinct") => ByRow(!ismissing))
end
island_extinct_names = map(get_species_names, island_extinct)
aggfactor = 8

rast = Raster(rand(10, 10, 5), (X(), Y(), Band([:a, :b, :c, :d, :e])))
stack = RasterStack(rast; layersfrom=Band)

lc_predictions_paths = (
    mus="$outputdir/lc_predictions_mus.nc",
    reu="$outputdir/lc_predictions_reu.nc",
    rod="$outputdir/lc_predictions_rod.nc",
)

# netcdf has the annoying center locus for time
lc_predictions = map(lc_predictions_paths) do path
    lc_predictions = RasterStack(path) |>
        x -> rebuild(Rasters.modify(BitArray, x); missingval=false) |>
        x -> Rasters.set(x, Ti => Int.(maybeshiftlocus(Start(), dims(x, Ti), )))
end

auxs = agg_aux(masks, slope_stacks, dems, lc_predictions, aggfactor)
Makie.image(lc_all[k][Ti=300]; colormap=:batlow)

# Just for Makie...
# This kind of merges the pixel colors...
lc_all = map(auxs) do aux
    rebuild(map(aux.lc) do lcs
        sum(map(.*, ntuple(UInt8, length(lcs)), lcs))
    end; missingval=0.0)
end
# Rasters.rplot(lc_all.mus; colorrange=(1, 6))
# Rasters.rplot(lc_predictions.mus[Ti=At(2000)])
# Rasters.rplot(lc_ag[Ti=At(2000)])
# Makie.plot(sum.(lc_fractions)[Ti=52])

k = :reu
k = :rod
k = :mus
include("species_rules.jl")
(; ruleset, rules, pred_ruleset, endemic_ruleset, islands) = def_syms(
    pred_df, introductions_df, island_extinct_tables, auxs, aggfactor; 
    replicates=nothing, pred_pops_aux
);
(; output, max_output, pred_output, init, output_kw) = islands[k]
@time sim!(max_output, pred_ruleset; proc=CPUGPU(), printframe=true);
@time sim!(output, ruleset; proc=SingleCPU(), printframe=true);
# @time sim!(pred_output, pred_ruleset; proc=CPUGPU(), printframe=true);
maxpops = maximum(max_output)
mkoutput = mk_pred(init, pred_ruleset; maxpops, landcover=lc_all[k], output_kw...)
mkoutput = mk_endemic(init, endemic_ruleset; maxpops, pred_pop=pred_pops_aux[k], landcover=lc_all[k], output_kw...)
mkoutput = mk(init, ruleset; maxpops, landcover=lc_all[k], output_kw...)
DynamicGrids.mask(mkoutput)
display(mkoutput)

# using ProfileView
# @profview sim!(output, ruleset; tspan=1550:1560, proc=SingleCPU(), printframe=true);
# using Cthulhu
# descend_clicked()
# sd = DynamicGrids.SimData(output, Ruleset(ruleset.rules[3]))
# @descend DynamicGrids.descendable(sd)
# Makie.plot(mean(map(i -> getproperty.(output[end-i].pred_pop, :macaque), 1:20)); colormap=:magma)
# maximum(getproperty.(mkoutput[end].pred_pops, :mouse))

# Get the max color for Makie
# max_pops = map(output) do frame
#     map(propertynames(first(frame.pred_pop))) do key
#         maximum(x -> getproperty(x, key), frame.pred_pop)
#     end
# end |> maximum

pred_pops_aux = map(islands) do island
    (; pred_output, init) = island
    @time sim!(pred_output, pred_ruleset; proc=SingleCPU(), printframe=true);
    A = cat(pred_output...; dims=3)
    DimArray(A, (dims(init.pred_pop)..., dims(pred_output)...))
end
pred_pops_aux.mus

function predict_extinctions(ruleset, islands)
    map(islands) do island
        _predict_extinctions(ruleset, island)
    end
end
function _predict_extinctions(ruleset, island)
    output = TransformedOutput(island.endemic_init; island.output_kw...) do f
        mean(sum, eachslice(f.endemic_presence; dims=3))
    end
    sim!(output, ruleset; printframe=true, proc=CPUGPU())
    return output
    # Return NamedVector of years to extinction
    # return sum(output) .+ first(island.output_kw.tspan) .- 1
end

(isnothing(pred_pops_aux) && error()); (; ruleset, pred_ruleset, endemic_ruleset, islands) = def_syms(
    pred_df, introductions_df, island_extinct_tables, auxs, aggfactor; 
    replicates=100, pred_pops_aux
);
(; output, max_output, pred_output, init, output_kw) = islands[k]
@time preds = predict_extinctions(endemic_ruleset, islands)
# @time sim!(output, endemic_ruleset; proc=CPUGPU(), printframe=true);
# @time sim!(max_output, ruleset; proc=SingleCPU(), printframe=true);
sum(map(xs -> xs .> 10, preds.reu)) .+ first(DynamicGrids.tspan(output))
sum(map(xs -> xs .> 10, preds.mus)) .+ first(DynamicGrids.tspan(output))
sum(map(xs -> xs .> 10, preds.rod)) .+ first(DynamicGrids.tspan(output))
