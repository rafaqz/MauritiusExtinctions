
using LossFunctions
using Optimization
using OptimizationOptimJL
using ThreadsX
using LatinHypercubeSampling

include("species_common.jl")
include("species_rules.jl")
include("makie.jl")

# Generate fakes

# function generate_fakes(island_keys, nspecies, nislandcounts; extant_extension=200)
#     fake_keys = Symbol.(Ref("sp"), 1:nspecies)
#     n = length(fake_keys)
#     Group = rand(["mammal", "bird", "reptile"], n)
#     # Mammals don't nest
#     Ground_nesting = map(Group) do g
#         g == "mammal" ? false : rand(Bool)
#     end
#     # Reptiles don't fly
#     Flight_capacity = map(Group) do g
#         g == "bird" ? rand(Bool) : true
#     end
#     fake_table = DataFrame(; Keys=fake_keys, Group, Ground_nesting, Flight_capacity)

#     no_replace_keys = copy(fake_keys)
#     island_fake_species = map(k -> Vector{Symbol}(undef, nislandcounts[k][k]), island_keys)
#     foreach(island_keys, island_fake_species) do k, sp
#         sample!(no_replace_keys, sp; replace=false)
#         setdiff!(no_replace_keys, sp)
#     end
#     onall = sample(no_replace_keys, nislandcounts[1].onall; replace=false)
#     setdiff!(no_replace_keys, onall)
#     mus_reu = sample(no_replace_keys, nislandcounts.mus.reu; replace=false)
#     setdiff!(no_replace_keys, mus_reu)
#     mus_rod = sample(no_replace_keys, nislandcounts.mus.rod; replace=false)
#     setdiff!(no_replace_keys, mus_rod)
#     reu_rod = sample(no_replace_keys, nislandcounts.reu.rod; replace=false)
#     setdiff!(no_replace_keys, reu_rod)
#     # Make sure we have used up all the keys
#     @assert isempty(no_replace_keys)
#     # Now add all the species endemic to more than one island to their islands
#     append!(island_fake_species.mus, mus_reu, mus_rod, onall)
#     append!(island_fake_species.reu, mus_reu, reu_rod, onall)
#     append!(island_fake_species.rod, mus_rod, reu_rod, onall)
#     # Make sure fake species exactly match pattern of endemism on real islands
#     @assert map(length, island_fake_species) == nislandspecies
#     island_fake_tables = map(island_fake_species) do keys
#         filter(r -> r.Keys in keys, fake_table)
#     end
#     EndemicNVs = map(island_fake_species) do fk
#         NamedVector{Tuple(fk),length(fk)}
#     end
#     (; fake_keys, island_fake_species, island_fake_tables, EndemicNVs)
# end

# Define the response parameters
# pred_keys = (:cat, :black_rat, :norway_rat, :mouse, :pig, :wolf_snake, :macaque)
pred_keys = (:cat, :black_rat, :pig, :norway_rat)
parameters = Model(predator_response_params(pred_keys))
range_times = [1970]
k = :mus
nreplicates = 48
# (; fake_keys, island_fake_species, island_fake_tables, EndemicNVs) = generate_fakes(island_keys, nspecies, nislandcounts; extant_extension)
# Load predator data
f = jldopen("sym_setup2_$aggfactor.jld2", "r")
auxs = f["auxs"];
pred_pops_aux = map(f["pred_pops_aux"]) do pp
    map(p -> p[pred_keys], pp)
end;
close(f)
pred_mean_density = mean(map(mean, pred_pops_aux)) |> NV{pred_keys,length(pred_keys)}
pred_max_density = map(1:length(first(first(pred_pops_aux)))) do i
    map(pred_pops_aux) do pred_pops
        maximum(getindex.(pred_pops, i))
    end |> mean
end |> NV{pred_keys,length(pred_keys)}
pred_carrycap = vcat(map(_ -> pred_max_density, parent(parameters))...)

killfactor = 0.65
randparams(x::AbstractArray{T}) where T = T.(rand(length(x)) .* killfactor)
# Assign random parameter values 0.0 .. 0.25, with low values more common
x0 = randparams(collect(parameters[:val]))
scale_params(x, cc) = x ./ sqrt.(cc) .* minimum(sqrt, cc)
# Scale for carrycap
parameters[:val] = scale_params(x0, pred_carrycap)
pred_response = parent(parameters)
pred_response

(; rules, endemic_ruleset, islands) = def_syms(
    pred_df, introductions_df, island_endemic_tables, auxs, aggfactor;
    replicates=nreplicates, pred_pops_aux, first_year, last_year, extant_extension,
    pred_keys, pred_response,
);
(; output, endemic_output, pred_output, init, output_kw) = islands[k];

# Calculate extinction dates for these parameters
p = (; endemic_ruleset, islands, last_year, pred_carrycap, extant_extension, loss=HuberLoss(), parameters, range_times);
endemics = extinction_forward(x0, p);
island_extinction_dates = map(e -> e.dates.mean, endemics);
map(minimum, island_extinction_dates)

# Calculate extinction ranges for these parameters
(; rules, endemic_ruleset, islands) = def_syms(
    pred_df, introductions_df, island_endemic_tables, auxs, aggfactor;
    replicates=nreplicates, pred_pops_aux, first_year, last_year, extant_extension,
    pred_keys, pred_response, island_extinction_dates,
);

island_ranges = map(islands) do island
    output = ArrayOutput(island.init; island.output_kw..., replicates=nothing)
    sim!(output, endemic_ruleset)
    return output[At([1970])]
end;


# Optimization
(; rules, endemic_ruleset, islands) = def_syms(
    pred_df, introductions_df, island_endemic_tables, auxs, aggfactor;
    replicates=nreplicates, pred_pops_aux, first_year, last_year, extant_extension,
    pred_keys, pred_response, island_extinction_dates,
);

# Save the original parameters
x_original = x0

# Plot the know ranges and make observables for predicted ranges
fig = Makie.Figure()
plot_obs = map(enumerate(island_ranges), island_extinction_dates) do (i, ranges), dates
    ax = Axis(fig[1, i]; ylabel="Target")
    obs_ax = Axis(fig[2, i]; ylabel="Predicted")
    dates_ax = Axis(fig[3, i]; ylabel="Extinction year", xzoomlock=true)
    dates = collect(dates)
    dates_obs = Observable(dates)
    endemics = ranges[1].endemic_presence
    diversity = Float64.(sum.(endemics))
    diversity_obs = Observable(diversity)
    colorrange = 0, length(first(endemics))
    image!(ax, diversity; colormap=:viridis, interpolate=false, colorrange)
    image!(obs_ax, diversity_obs; colormap=:viridis, interpolate=false, colorrange)
    bp1 = Makie.barplot!(dates_ax, dates)
    bp2 = Makie.barplot!(dates_ax, dates_obs)
    Makie.ylims!(dates_ax; low=1500, high=2200)
    return (; diversity_obs, dates_obs)
end |> NamedTuple{keys(islands)}
loss_obs = Observable("loss")
params_obs = Observable(x_original)
params_ax = Axis(fig[4, 1:3];
    xticks=axes(collect(parameters[:val]), 1),
    xtickformat=I -> map(string, collect(parameters[:label])[map(Int, I)]),
    xticklabelrotation=45.0,
    xlabel=loss_obs,
    ylabel="Parameter val",
    xzoomlock=true,
)
Makie.barplot!(params_ax, x_original)
Makie.barplot!(params_ax, params_obs)

# Define new random parameters
x0 = randparams(collect(parameters[:val]))
lb = map(x0 -> zero(x0), x0)
ub = map(x0 -> oneunit(x0), x0)
throwit = Ref(false)
obs = (; plot_obs, params_obs, loss_obs, throwit)

p = (; endemic_ruleset, islands, parameters, last_year, extant_extension, loss=HuberLoss(), range_times, island_ranges, pred_carrycap, obs);
extinction_objective(x0, p);
# Burn in std

# Optimization.jl
prob = OptimizationProblem(extinction_objective, x0, p; lb, ub)
sol = solve(prob, SAMIN(; verbosity=3, rt=0.8); maxiters=20000)
# sol = solve(prob, ParticleSwarm(lower=lb, upper=ub, n_particles=25); maxiters=2500)

# Optim.jl
# optim_f = let p = p
#     x -> extinction_objective(x, p)
# end
# res = Optim.optimize(optim_f, lb, ub, x0, SAMIN(; verbosity=3))

# Analysis
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

endemic_ruleset[:val] = scale_params(x0, pred_carrycap)
(; rules, endemic_ruleset, islands) = def_syms(
    pred_df, introductions_df, island_endemic_tables, auxs, aggfactor;
    replicates=nothing, pred_pops_aux, first_year, last_year, extant_extension,
    pred_keys, island_extinction_dates, 
);
k = :mus
(; output, endemic_output, pred_output, init, output_kw) = islands[k];
tspan = first_year:last_year
mk_endemic(init, endemic_ruleset; landcover=lc_all[k], pred_pop=pred_pops_aux[k], tspan, islands[k].output_kw...)


# Hyper cube sampling. Plan optimization takes a long time.
nsimulations = 1_000
nparameters = 20#length(parameters[:val])
ngenerations = 10
plan, fitness = LHCoptim(nsimulations, nparameters, ngenerations)#; threading=true)
plan
fitness
scaled_plan = scaleLHC(plan, [(0.0, 1.0) for _ in 1:20])

for i in 1:nsimulations
    x = plan[100, :]
    loss = extinction_forward(x, p)
end
