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
using Setfield
using ConstructionBase
using JLD2

# includet("optimisation.jl")
include("species_rules.jl")
include("species_tables.jl")
include("raster_common.jl")
include("makie.jl")

uncleared = gpu_cleanup(Rasters.modify(BitArray, Raster("uncleared.nc")))
# forested = gpu_cleanup(modify(BitArray, Raster("forested.nc")))
# Set aggregation
aggfactor = 8
# And the last year of the simulation
last_year = 2018

# Build auxiliary rasters
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
Makie.plot(masks.rod)
# Remove islands of rodrigues
masks.rod[X=60 .. 63.33, Y=-19.775 .. -19.675] .= false
masks.rod[X=63 .. 65, Y= -19.8 .. -19.775] .= false
auxs = agg_aux(masks, slope_stacks, dems, lc_predictions, aggfactor, last_year)
# Just for Makie, kind of merges the pixel colors...
lc_all = map(auxs) do aux
    rebuild(map(aux.lc) do lcs
        sum(map(.*, ntuple(UInt8, length(lcs)), lcs))
    end; missingval=0.0)
end

# Set up and run simulations
k = :reu
k = :rod
k = :mus
include("species_rules.jl")
pred_pops_aux = isdefined(Main, :pred_pops_aux) ? pred_pops_aux : map(_ -> nothing, dems)
(; ruleset, rules, pred_ruleset, endemic_ruleset, islands, pred_response) = def_syms(
    pred_df, introductions_df, island_endemic_tables, auxs, aggfactor; 
    replicates=nothing, pred_pops_aux
);
(; output, endemic_output, pred_output, init, output_kw) = islands[k]
# @time sim!(output, ruleset; proc=SingleCPU(), printframe=true);
# Optimisation
pred_pops_aux = map(islands) do island
    (; pred_output, init) = island
    @time sim!(pred_output, pred_ruleset; proc=SingleCPU(), printframe=true);
    A = cat(pred_output...; dims=3)
    DimArray(A, (dims(init.pred_pop)..., dims(pred_output)...))
end
(; ruleset, rules, pred_ruleset, endemic_ruleset, islands, pred_response) = def_syms(
    pred_df, introductions_df, island_endemic_tables, auxs, aggfactor; 
    replicates=2, pred_pops_aux
);



# Store so we don't have to run the above
jldsave("sym_setup.jld2";
    pred_df, introductions_df, island_endemic_tables, auxs, aggfactor, pred_pops_aux, 
    ruleset, rules, pred_ruleset, endemic_ruleset, islands, pred_response
);

(; ruleset, rules, pred_ruleset, endemic_ruleset, islands) = def_syms(
    pred_df, introductions_df, island_endemic_tables, auxs, aggfactor; 
    replicates=nothing, pred_pops_aux
);
(; output, endemic_output, pred_output, init, output_kw) = islands[k]

@time sim!(endemic_output, endemic_ruleset; proc=SingleCPU(), printframe=true);

mkoutput = mk(init, ruleset; landcover=lc_all[k], output_kw..., ncolumns=5)
mkoutput = mk_pred(init, pred_ruleset; landcover=lc_all[k], output_kw...)
sim!(mkoutput, ruleset; proc=SingleCPU(), printframe=true);
display(mkoutput)

# k = :mus
# k = :rod
# k = :reu
# (; output, endemic_output, pred_output, init, output_kw) = islands[k]

# p = Rasters.rplot(lc_all.mus[Ti=At(1700:2018)]; colorrange=(1, 6))
# save("images/landcover_simulation.png", p)
# # mkoutput = mk_endemic(init, endemic_ruleset; ncolumns=5, maxpops, pred_pop=pred_pops_aux[k], landcover=lc_all[k], output_kw...)
# display(mkoutput)


@time sim!(endemic_output, endemic_ruleset; proc=SingleCPU(), printframe=true);
cu_endemic_output = Adapt.adapt(CuArray, endemic_output)
CUDA.@profile sim!(cu_endemic_output, endemic_ruleset; proc=CuGPU(), printframe=true)
@time preds = predict_extinctions(endemic_ruleset, islands)
# @time sim!(output, endemic_ruleset; proc=CPUGPU(), printframe=true);
sum(map(xs -> xs .> 10, preds.reu)) .+ first(DynamicGrids.tspan(output))

# intercept = 0.085
# slope = 1.177
# predmass = 1000
# preymass = 1000 * slope + intercept
1:0.1:log(1000)

using Distributions
n = Distributions.Normal(41.0, 102.0)
Makie.plot(n)
scalar = 1 / Distributions.pdf(n, 41)
Distributions.pdf(n, 16) * scalar
Distributions.pdf(n, 220) * scalar
Distributions.pdf(n, 340) * scalar

