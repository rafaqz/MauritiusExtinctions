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
@async run(`libreoffice $mascarene_species_csv`)
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
@time sim!(max_output, pred_ruleset; proc=SingleCPU(), printframe=true);
# @time sim!(output, ruleset; proc=SingleCPU(), printframe=true);
# Optimisation
pred_pops_aux = map(islands) do island
    (; pred_output, init) = island
    @time sim!(pred_output, pred_ruleset; proc=SingleCPU(), printframe=true);
    A = cat(pred_output...; dims=3)
    DimArray(A, (dims(init.pred_pop)..., dims(pred_output)...))
end
(; ruleset, rules, pred_ruleset, endemic_ruleset, islands) = def_syms(
    pred_df, introductions_df, island_endemic_tables, auxs, aggfactor; 
    replicates=nothing, pred_pops_aux
);
(; output, max_output, endemic_output, pred_output, init, output_kw) = islands[k]
island = islands[k]

endemic_ruleset = Ruleset(
    Chain(rules.endemic_recouperation_rule, rules.aux_pred_risks_rule, rules.clearing_rule);
    boundary=Remove()
)
@time sim!(endemic_output, endemic_ruleset; proc=SingleCPU(), printframe=true);
# using ProfileView
# using CUDA, Adapt
# cu_endemic_output = Adapt.adapt(CuArray, endemic_output)
# CUDA.@profile sim!(cu_endemic_output, endemic_ruleset; proc=CuGPU(), printframe=true)
# DynamicGrids.aux(cu_endemic_output).pred_effect

# @time sim!(output, endemic_ruleset; proc=CuGPU(), printframe=true);
# @profview sim!(output, endemic_ruleset; proc=SingleCPU(), printframe=true);
maxpops = maximum(max_output)
maxpops = maximum(pred_pops_aux[k])

mkoutput = mk_pred(init, pred_ruleset; maxpops, landcover=lc_all[k], output_kw...)
mkoutput = mk_endemic(init, endemic_ruleset; ncolumns=5, maxpops, pred_pop=pred_pops_aux[k], landcover=lc_all[k], output_kw...)
mkoutput = mk(init, ruleset; maxpops, landcover=lc_all[k], output_kw..., ncolumns=5)
display(mkoutput)

function predict_extinctions(ruleset, islands, parameters)
    map(islands) do island
        _predict_extinctions(ruleset, island, parameters)
    end
end
function _predict_extinctions(ruleset, island, pred_pops, pred_response, parameters)
    aux = island.output_kw.aux
    traits = aux.endemic_traits
    pred_suscept = predator_suceptibility(island.mass_response, pred_response, traits)
    pred_pop = pred_pops_aux[k];
    pred_pop
    using ProfileView, Cthulhu
    @descend 
    f(parent(pred_pop), pred_suscept);
    f(pred_pop, pred_suscept) = 
    predator_effect.(tanh, pred_pop, (pred_suscept,))
    typeof(parent(pred_pop))
    aux = (;
        recouperation_rate=island.recouperation_rates,
        habitat_requirement=island.habitat_requirements,
        pred_pop,
        pred_effect,
    )

    output = TransformedOutput(island.endemic_init; island.output_kw..., aux) do f
        mean(sum, eachslice(f.endemic_presence; dims=3))
    end
    sim!(output, ruleset; printframe=true, proc=CPUGPU())
    return output
    # Return NamedVector of years to extinction
    # return sum(output) .+ first(island.output_kw.tspan) .- 1
end

function predator_effect(f, pred_pop, pred_suscept)
    Float32.(f.(sum(map(pred_suscept, pred_pop) do ps, pp
        map(xs -> xs .* pp, ps)
    end)))
end

function predator_suceptibility(mass_response, pred_response, traits)
    pred_suscept = mapreduce(+, pred_response, traits) do pr, t
        map(mass_response, pr) do m, p
            t .* p .* m
        end
    end ./ (32 * 8^2)
end

# Outputs with replicates
(isnothing(pred_pops_aux) && error()); (; ruleset, pred_ruleset, endemic_ruleset, islands) = def_syms(
    pred_df, introductions_df, island_endemic_tables, auxs, aggfactor; 
    replicates=10, pred_pops_aux
);
(; output, max_output, endemic_output, pred_output, init, output_kw) = islands[k]
@time sim!(endemic_output, endemic_ruleset; proc=SingleCPU(), printframe=true);
cu_endemic_output = Adapt.adapt(CuArray, endemic_output)
CUDA.@profile sim!(cu_endemic_output, endemic_ruleset; proc=CuGPU(), printframe=true)
@time preds = predict_extinctions(endemic_ruleset, islands)
# @time sim!(output, endemic_ruleset; proc=CPUGPU(), printframe=true);
# @time sim!(max_output, ruleset; proc=SingleCPU(), printframe=true);
sum(map(xs -> xs .> 10, preds.reu)) .+ first(DynamicGrids.tspan(output))
sum(map(xs -> xs .> 10, preds.mus)) .+ first(DynamicGrids.tspan(output))
sum(map(xs -> xs .> 10, preds.rod)) .+ first(DynamicGrids.tspan(output))

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

