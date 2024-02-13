using ModelParameters
using StaticArrays
using DynamicGrids
using Dispersal
using LandscapeChange
using DataFrames
using TerminalPager
using Rasters
using Statistics
using Setfield
using ConstructionBase
using JLD2
using LossFunctions
using Optimization
using OptimizationOptimJL

function predict_extinctions(endemic_ruleset::Ruleset, islands, pred_response)
    map(islands) do island
        _predict_extinctions(endemic_ruleset, island, pred_response)
    end
end
function _predict_extinctions(endemic_ruleset::Ruleset, island, pred_response)
    kw = island.output_kw
    aux = kw.aux
    pred_suscept = predator_suceptibility(pred_response, aux.endemic_traits)
    (; pred_pop, pred_effect) = aux
    generate_predator_effect!(tanh, pred_effect, pred_pop, pred_suscept)
    output = TransformedOutput(island.endemic_init; kw...) do f
        # Take the mean over the replicates dimension
        mean(sum, eachslice(f.endemic_presence; dims=3))
    end
    sim!(output, endemic_ruleset; printframe=true, proc=CPUGPU())
    return output
end

function extinction_objective(x, p)
    (; endemic_ruleset, islands, parameters, loss) = p
    # Update
    parameters[:val] = x
    pred_response = stripparams(parameters)
    predictions = predict_extinctions(endemic_ruleset, islands, pred_response)
    island_losses = map(predictions, islands) do preds, island
        firstyear = first(island.output_kw.tspan)
        years_present_per_species = sum(preds)
        losses = loss.(island.extinction_dates, years_present_per_species .+ firstyear)
        q = quantile(losses, 0.9)
        @show q
        sum(losses) do l
            l > q ? zero(l) : l
        end
    end
    return sum(island_losses)
end

function extinction_forward(x, p)
    (; endemic_ruleset, islands, parameters, loss) = p
    # Update
    parameters[:val] = x
    pred_response = stripparams(parameters)
    predict_extinctions(endemic_ruleset, islands, pred_response)
end

f = jldopen("sym_setup.jld2", "r") 
pred_df = f["pred_df"] 
introductions_df = f["introductions_df"] 
island_endemic_tables = f["island_endemic_tables"] 
aggfactor = f["aggfactor"] 
pred_pops_aux = f["pred_pops_aux"] 
auxs = f["auxs"] 
endemic_ruleset = f["endemic_ruleset"]
islands = f["islands"]
pred_response = f["pred_response"]

# endemic_ruleset = Ruleset(
#     Chain(rules.endemic_recouperation_rule, rules.aux_pred_risks_rule, rules.clearing_rule);
#     boundary=Remove()
# )
p = (; 
    endemic_ruleset, 
    islands, 
    loss=HuberLoss(),
    parameters=Model(pred_response)
);
x0 = collect(p.parameters[:val])
predictions = extinction_objective(x0, p)

# Outputs with replicates
include("species_rules.jl")
(isnothing(pred_pops_aux) && error()); (; endemic_ruleset, islands) = def_syms(
    pred_df, introductions_df, island_endemic_tables, auxs, aggfactor; 
    replicates=10, pred_pops_aux
);
(; output, max_output, endemic_output, pred_output, init, output_kw) = islands[k]
predictions = predict_extinctions(endemic_ruleset, islands, pred_response)
keys(predictionsj)

prob = OptimizationProblem(extinction_objective, x0, p)
sol = solve(prob, NelderMead())

641629.6, = 1.11283155e6
q = 31250.5

q = 336834.5000000041
q = 1.6668971500000039e6
q = 13356.1

q = 443735.4
q = 991893.9000000001
q = 24152.9

p.parameters
x = sol.u
preds = extinction_forward(x, p)
preds.reu[At(1600)] |> NamedTuple |> pairs

pred_pop = pred_pops_aux.mus
maxpops = maximum(pred_pop)
tspan = islands.mus.output_kw.tspan
kw = islands.mus.output_kw
landcover = lc_all.mus

(isnothing(pred_pops_aux) && error()); (; endemic_ruleset, islands) = def_syms(
    pred_df, introductions_df, island_endemic_tables, auxs, aggfactor; 
    replicates=nothing, pred_pops_aux
);
mk_endemic(init, ruleset; maxpops, landcover, pred_pop, tspan, kw...)
