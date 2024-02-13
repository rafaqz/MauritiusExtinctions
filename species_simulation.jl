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
include("raster_common.jl")
include("makie.jl")

uncleared = gpu_cleanup(Rasters.modify(BitArray, Raster("uncleared.nc")))
# forested = gpu_cleanup(modify(BitArray, Raster("forested.nc")))

# Import tabular data
pred_df = CSV.read("tables/animals.csv", DataFrame)
introductions_df = CSV.read("tables/introductions.csv", DataFrame)
mascarene_species_csv = "tables/mascarene_species.csv"
# @async run(`libreoffice $mascarene_species_csv`)
all_species = CSV.read(mascarene_species_csv, DataFrame) |> 
    x -> subset(x, :Species => ByRow(!ismissing))
island_tables = map(island_keys) do key
    df = DataFrame(subset(all_species, key => x -> .!ismissing.(x)))
    df.extinct = map(df[!, "$(key)_extinct"]) do e
        ismissing(e) ? missing : eval(Meta.parse(e))::UnitRange
    end
    df.introduced = map(df[!, "$(key)_introduced"]) do e
        ismissing(e) ? missing : eval(Meta.parse(e))::UnitRange
    end
    df
end
island_endemic_tables = map(island_tables) do tbl
    # TODO add missing mass rows and remove the missing Mass filter
    DataFrame(subset(tbl, :Origin => ByRow(==("Endemic")), :Mass => ByRow(!ismissing); skipmissing=true))
end
get_species_names(table) = Tuple(Symbol.(replace.(skipmissing(table.Species), Ref(" " => "_"))))
island_names = NamedTuple{keys(island_tables)}(keys(island_tables))
island_endemic_names = map(island_tables, island_names) do table, name
    get_species_names(subset(table, Symbol(name, "_extinct") => ByRow(!ismissing)))
end

# Set aggregation
aggfactor = 8

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
auxs = agg_aux(masks, slope_stacks, dems, lc_predictions, aggfactor)
auxs.mus.lc
lc = map(CartesianIndices(first(lc_ag1))) do I
    Float32.(NV(lc_ag1[I]))
end
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
(; output, max_output, endemic_output, pred_output, init, output_kw) = islands[k]
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

# (; output, max_output, endemic_output, pred_output, init, output_kw) = islands[k]
# @time sim!(max_output, pred_ruleset; proc=SingleCPU(), printframe=true);
# # maxpops = maximum(max_output)

# (; ruleset, rules, pred_ruleset, endemic_ruleset, islands) = def_syms(
#     pred_df, introductions_df, island_endemic_tables, auxs, aggfactor; 
#     replicates=nothing, pred_pops_aux
# );
# (; output, max_output, endemic_output, pred_output, init, output_kw) = islands[k]
# island = islands[k]

# @time sim!(endemic_output, endemic_ruleset; proc=SingleCPU(), printframe=true);

# mkoutput = mk_pred(init, pred_ruleset; maxpops, landcover=lc_all[k], output_kw...)
# mkoutput = mk(init, ruleset; maxpops, landcover=lc_all[k], output_kw..., ncolumns=5)

# k = :mus
# k = :rod
# k = :reu
# (; output, max_output, endemic_output, pred_output, init, output_kw) = islands[k]

# p = Rasters.rplot(lc_all.mus[Ti=At(1700:2018)]; colorrange=(1, 6))
# save("images/landcover_simulation.png", p)
# # mkoutput = mk_endemic(init, endemic_ruleset; ncolumns=5, maxpops, pred_pop=pred_pops_aux[k], landcover=lc_all[k], output_kw...)
# display(mkoutput)


@time sim!(endemic_output, endemic_ruleset; proc=SingleCPU(), printframe=true);
cu_endemic_output = Adapt.adapt(CuArray, endemic_output)
CUDA.@profile sim!(cu_endemic_output, endemic_ruleset; proc=CuGPU(), printframe=true)
@time preds = predict_extinctions(endemic_ruleset, islands)
# @time sim!(output, endemic_ruleset; proc=CPUGPU(), printframe=true);
# @time sim!(max_output, ruleset; proc=SingleCPU(), printframe=true);
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

