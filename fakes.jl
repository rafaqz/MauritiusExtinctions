include("species_common.jl")

using LossFunctions
using Optimization
using OptimizationOptimJL
using ThreadsX

# Load predator data

f = jldopen("sym_setup.jld2", "r")
pred_pops_aux = f["pred_pops_aux"];
pred_response = f["pred_response"];
auxs = f["auxs"];

# Generate fakes

function generate_fakes(island_keys, nspecies, nislandcounts; extant_extension=200)
    fake_keys = Symbol.(Ref("sp"), 1:nspecies)
    n = length(fake_keys)
    Group = rand(["mammal", "bird", "reptile"], n)
    Ground_nesting = rand(Bool, n)
    Flight_capacity = rand(Bool, n)
    fake_table = DataFrame(; Keys=fake_keys, Group, Ground_nesting, Flight_capacity) 

    no_replace_keys = copy(fake_keys)
    island_fake_species = map(k -> Vector{Symbol}(undef, nislandcounts[k][k]), island_keys) 
    foreach(island_keys, island_fake_species) do k, sp
        sample!(no_replace_keys, sp; replace=false)
        setdiff!(no_replace_keys, sp) 
    end
    onall = sample(no_replace_keys, nislandcounts[1].onall; replace=false)
    setdiff!(no_replace_keys, onall) 
    mus_reu = sample(no_replace_keys, nislandcounts.mus.reu; replace=false)
    setdiff!(no_replace_keys, mus_reu) 
    mus_rod = sample(no_replace_keys, nislandcounts.mus.rod; replace=false)
    setdiff!(no_replace_keys, mus_rod) 
    reu_rod = sample(no_replace_keys, nislandcounts.reu.rod; replace=false)
    setdiff!(no_replace_keys, reu_rod) 
    # Make sure we have used up all the keys
    @assert isempty(no_replace_keys)
    # Now add all the species endemic to more than one island to their islands
    append!(island_fake_species.mus, mus_reu, mus_rod, onall)
    append!(island_fake_species.reu, mus_reu, reu_rod, onall)
    append!(island_fake_species.rod, mus_rod, reu_rod, onall)
    # Make sure fake species exactly match pattern of endemism on real islands
    @assert map(length, island_fake_species) == nislandspecies
    island_fake_tables = map(island_fake_species) do keys
        filter(r -> r.Keys in keys, fake_table)
    end
    EndemicNVs = map(island_fake_species) do fk
        NamedVector{Tuple(fk),length(fk)}
    end
    (; fake_keys, island_fake_species, island_fake_tables, EndemicNVs)
end

(; fake_keys, island_fake_species, island_fake_tables, EndemicNVs) = generate_fakes(island_keys, nspecies, nislandcounts; extant_extension)

# Define the response parameters
pred_keys = (:cat, :black_rat, :norway_rat, :mouse, :pig, :wolf_snake, :macaque)

k = :mus
include("species_rules.jl")
nreplicates = 25
island_extinction_dates = isdefined(Main, :island_extinction_dates) ? island_extinction_dates : map(_ -> nothing, EndemicNVs)
(isnothing(pred_pops_aux) && error()); (; rules, endemic_ruleset, islands) = def_syms(
    pred_df, introductions_df, island_fake_tables, auxs, aggfactor;
    replicates=nreplicates, pred_pops_aux, last_year, extant_extension,
    pred_keys, pred_response, EndemicNVs, # Other fake bits
    island_extinction_dates, island_mass_response = map(_ -> nothing, EndemicNVs),
);
(; output, endemic_output, pred_output, init, output_kw) = islands[k];


# Calculate extinction dates for these parameters

parameters = Model(predator_response_params(pred_keys))
# Assign random parameter values 0.0 .. 0.25, with low values more common
randparams!(parameters) = parameters[:val] = Float32.((rand(length(parameters[:val])) ./ 5) .^ 2)
randparams!(parameters)
pred_response = parent(parameters)
# We don't know the dates yet
x0 = collect(parameters[:val])
lb = map(x0 -> zero(x0), x0)
ub = map(x0 -> oneunit(x0), x0)

p = (; endemic_ruleset, islands, last_year, extant_extension, loss=HuberLoss(), parameters);
endemics = extinction_forward(x0, p)
island_extinction_dates = map(e -> e.dates.mean, endemics)
# tspan = islands.mus.output_kw.tspan
# species = collect(propertynames(first(endemics.mus.dates.years)))
# m = DimMatrix(reduce(hcat, endemics.mus.dates.years), (; species, replicates=1:nreplicates))


# Optimization 

(isnothing(pred_pops_aux) && error()); (; rules, endemic_ruleset, islands) = def_syms(
    pred_df, introductions_df, island_fake_tables, auxs, aggfactor;
    replicates=nreplicates, pred_pops_aux, last_year, extant_extension,
    pred_keys, pred_response, EndemicNVs, # Other fake bits
    island_extinction_dates, island_mass_response = map(_ -> nothing, EndemicNVs),
);
# New random parameters
x_original = collect(p.parameters[:val])
randparams!(p.parameters)
x0 = collect(p.parameters[:val])
p = (; endemic_ruleset, islands, last_year, extant_extension, loss=HuberLoss(), parameters);
@time extinction_objective(x0, p)
ub = fill(1.0, 35)
lb = fill(0.0, 35)
x0 = rand(35)
prob = OptimizationProblem(extinction_objective, x0, p; lb, ub)
@time sol = solve(prob, NelderMead())

using Pkg; Pkg.status(; mode = PKGMODE_MANIFEST)

# Analisys

parameters[:val] = x_original
p = (; endemic_ruleset, islands, last_year, extant_extension, loss=HuberLoss(), parameters);
original_run = extinction_forward(x0, p)
original_extinction_dates = map(e -> e.dates.mean, original_run)

parameters[:val] = sol.u
p = (; endemic_ruleset, islands, last_year, extant_extension, loss=HuberLoss(), parameters);
optimised_run = extinction_forward(x0, p)
optimised_extinction_dates = map(e -> e.dates.mean, optimised_run)

param_comparison = x_original .=> sol.u

NamedTuple{p.parameters[:label]}((x_original)) |> pairs
NamedTuple{p.parameters[:label]}((sol.u)) |> pairs
mean(map(abs, NamedTuple{p.parameters[:label]}((x_original .- sol.u))))
mean(map(abs, sol.u))
A |> pairs

original_extinction_dates.mus .=> optimised_extinction_dates.mus
original_extinction_dates.reu .=> optimised_extinction_dates.reu
original_extinction_dates.rod .=> optimised_extinction_dates.rod

original_extinction_dates.reu .- optimised_extinction_dates.reu
original_extinction_dates.rod .- optimised_extinction_dates.rod

diff = reduce(vcat, original_extinction_dates) .- reduce(vcat, optimised_extinction_dates)
mean(abs.(diff))

endemics.mus.dates.years
endemics.reu.dates.years
endemics.rod.dates.years
endemics.mus.dates.mean
endemics.reu.dates.mean
endemics.rod.dates.mean
endemics.mus.dates.std
endemics.reu.dates.std
endemics.rod.dates.std

mean(vcat(endemics.mus.dates.std, endemics.reu.dates.std, endemics.rod.dates.std))
m = DimMatrix(reduce(hcat, endemics.mus.dates.years), (; species=island_fake_keys.mus, replicates=1:nreplicates))

include("species_rules.jl")
(isnothing(pred_pops_aux) && error()); (; rules, endemic_ruleset, islands) = def_syms(
    pred_df, introductions_df, island_fake_tables, auxs, aggfactor;
    replicates=nothing, pred_pops_aux, last_year, extant_extension,
    pred_keys, pred_response, EndemicNVs, # Other fake bits
    island_extinction_dates, island_mass_response = map(_ -> nothing, EndemicNVs),
);
(; output, endemic_output, pred_output, init, output_kw) = islands[k];
k = :mus
mk_endemic(init, endemic_ruleset; landcover=lc_all[k], pred_pop=pred_pops_aux[k], tspan, ncolumns=5, islands[k].output_kw...)
 

using Optim

rosenbrock(x, p) = (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
x0 = zeros(2)
lb = fill(-1.0, 2)
ub = fill(200.0, 2)
_p = [1.0, 100.0]

l1 = rosenbrock(x0, _p)
prob = OptimizationProblem(rosenbrock, x0, _p; lb, ub)
solve(prob, NelderMead())
