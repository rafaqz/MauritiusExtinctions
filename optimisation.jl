using ModelParameters
using StaticArrays
using DynamicGrids
using Dispersal
using LandscapeChange
using Stencils
using DataFrames
using TerminalPager
using Rasters
using Statistics
using Setfield
using ConstructionBase
using LossFunctions
using Optimization
using OptimizationOptimJL
using ThreadsX
using GLMakie
using NCDatasets
using Geomorphometry
using JLD2

f = jldopen("sym_setup.jld2", "r")
pred_df = f["pred_df"];
introductions_df = f["introductions_df"];
island_endemic_tables = f["island_endemic_tables"];
aggfactor = f["aggfactor"];
pred_pops_aux = f["pred_pops_aux"];
auxs = f["auxs"];
pred_response = f["pred_response"];

k = :mus
k = :reu
k = :rod
# Outputs with replicates
include("species_rules.jl")
(isnothing(pred_pops_aux) && error()); (; endemic_ruleset, islands) = def_syms(
    pred_df, introductions_df, island_endemic_tables, auxs, aggfactor;
    replicates=10, pred_pops_aux
);
(; output, endemic_output, pred_output, init, output_kw) = islands[k];

@time predictions = predict_timeline(endemic_ruleset, islands, pred_response);

using ProfileView
predictions = predict_timeline(endemic_ruleset, islands, pred_response);

p = (;
    endemic_ruleset,
    islands,
    loss=HuberLoss(),
    parameters=Model(pred_response)
);
x0 = collect(p.parameters[:val])
lb = first.(p.parameters[:bounds])
ub = last.(p.parameters[:bounds])
prob = OptimizationProblem(extinction_objective, x0, p; lb, up)
@time sol = solve(prob, NelderMead())
x = sol.u
endemics = extinction_forward(x, p)
endemics.mus.dates.years |> pairs


pred_pop = pred_pops_aux[k]
tspan = islands[k].output_kw.tspan
kw = islands[k].output_kw
landcover = lc_all[k]

p.parameters[:label] .=> p.parameters[:val] |> collect

# Run the result visually
(isnothing(pred_pops_aux) && error()); (; endemic_ruleset, islands) = def_syms(
    pred_df, introductions_df, island_endemic_tables, auxs, aggfactor;
    replicates=nothing, pred_pops_aux
);
(; output, max_output, endemic_output, pred_output, init, output_kw) = islands[k]
endemic_ruleset[:val] = sol.u
mk_endemic(init, endemic_ruleset; landcover, pred_pop, tspan, ncolumns=5, kw...)
