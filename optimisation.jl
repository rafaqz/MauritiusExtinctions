using LossFunctions
using Optimization
using OptimizationOptimJL
using ThreadsX
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

f = jldopen("sym_setup.jld2", "r")
pred_df = f["pred_df"];
introductions_df = copy(f["introductions_df"]);
aggfactor = f["aggfactor"];
pred_pops_aux = f["pred_pops_aux"];
auxs = f["auxs"];
pred_response = f["pred_response"];
# Run the result visually
(isnothing(pred_pops_aux) && error()); (; endemic_ruleset, islands) = def_syms(
    pred_df, introductions_df, island_endemic_tables, auxs, aggfactor;
    replicates=nothing, pred_pops_aux
);
(; output, endemic_output, pred_output, init, output_kw) = islands[k]

last_year = 2018
extant_extension = 200

k = :mus
k = :reu
k = :rod
# Outputs with replicates
include("species_rules.jl")
map(island_endemic_tables) do table
    sort!(table, :GBIFSpecies)
end;
(isnothing(pred_pops_aux) && error()); (; endemic_ruleset, islands) = def_syms(
    pred_df, introductions_df, island_endemic_tables, auxs, aggfactor;
    replicates=25, pred_pops_aux, last_year, extant_extension
);
(; output, endemic_output, pred_output, init, output_kw) = islands[k];

@time predictions = predict_timeline(endemic_ruleset, islands, pred_response);


# Optimization
p = (;
    endemic_ruleset,
    islands,
    last_year,
    extant_extension,
    loss=HuberLoss(),
    parameters=Model(pred_response)
);
x0 = collect(p.parameters[:val])
lb = map(x -> zero(x), x0)
ub = map(x -> oneunit(x), x0)

@time extinction_objective(x0, p)

prob = OptimizationProblem(extinction_objective, x0, p)#; lb, ub)
@time sol = solve(prob, NelderMead())

p.parameters[:label] .=> p.parameters[:val] |> collect

# Resultsendemics
x = sol.u
p.parameters[:val] = x
vals = collect(p.parameters[:val])
labels = collect(p.parameters[:label])
endemics = extinction_forward(x, p)
jldsave("sol"; sol);



e = jldopen("endemic.jld2", "r")
endemics = e["endemics"];

# Analysis
# plt = Makie.barplot(vals)
# Makie.xlabel!(plt.axis, string.(label))

results = map(map(deepcopy, island_endemic_tables), islands, endemics) do table, island, endmic
    sort!(table, :GBIFSpecies)
    keys = first.(sort(pairs(NamedTuple(island.extinction_dates))))
    extinction = last.(sort(pairs(NamedTuple(island.extinction_dates))))
    predicted = last.(sort(pairs(NamedTuple(endmic.dates.mean))))
    insertcols!(table, 1, :key=>keys, :target=>extinction, :predicted=>predicted)
    table
end
function plot_results(fig, i, res)
    axis = Axis(fig[1, i])
    sort!(res, :Mass)
    sort!(res, :Group)
    x = 1.0:length(res.target)
    y = mean.(zip(res.target, res.predicted))
    yerr = (res.target .- res.predicted) ./ 2
    errcolor = (x -> x > 0 ? :red : :blue).(yerr)
    p = errorbars!(axis, x, y, yerr; whiskerwidth = 12, color=errcolor)
    Makie.scatter!(axis, res.target; color=:black, marker=:xcross, markersize=20)
    Makie.scatter!(axis, parent(parent(res.predicted)); color=:red, markersize=20)
    text = string.(res.GBIFSpecies, :_predicted)
    Makie.text!(axis, res.predicted; rotation=0.6, text)
    text = string.(res.GBIFSpecies)
    Makie.text!(axis, res.target; rotation=-0.6, text)
end
fig = Figure()
plot_results.(Ref(fig), 1:3, collect(results))
names(island_endemic_tables.mus)


map(tuple, endemics[k].dates.mean, islands[k].extinction_dates) |> NamedTuple |> pairs |> sort
map(tuple, endemics[k].dates.mean, islands[k].extinction_dates) |> NamedTuple |> pairs |> sort
map(-, endemics[k].dates.mean, islands[k].extinction_dates) |> NamedTuple |> pairs |> sort

# Run the result visually
(isnothing(pred_pops_aux) && error()); (; endemic_ruleset, islands) = def_syms(
    pred_df, introductions_df, island_endemic_tables, auxs, aggfactor;
    replicates=nothing, pred_pops_aux
);
(; output, endemic_output, pred_output, init, output_kw) = islands[k]

pred_pop = pred_pops_aux[k]
tspan = islands[k].output_kw.tspan
kw = islands[k].output_kw
landcover = lc_all[k]
endemic_ruleset[:val] = sol.u
mk_endemic(init, endemic_ruleset; landcover, pred_pop, tspan, ncolumns=5, kw...)
