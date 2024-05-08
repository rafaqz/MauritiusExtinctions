using Graphs, GLMakie, GraphMakie, NetworkLayout, InvertedIndices, LandscapeChange
using CSV, DataFrames, StaticArrays

# convenience functions
# picks n random numbers between low and high, that are no more than a factor 2 apart
function rejection_range(n, low=0, high=1, max_multiplier=2)
    while true
        cont = false
        r = rand(n) .* (high - low) .+ low
        for i in 1:n, j in 1:n
            if r[i] / r[j] > max_multiplier || r[j] / r[i] > max_multiplier
                cont = true
            end
        end
        cont || return r
    end
end

# visualize a single interaction matrix
function plot_interactions(mat, nodes)
    vals = [mat[j, i] for i in axes(mat, 1), j in axes(mat, 2) if i != j]
    cols = ifelse.(vals .>= 0, :blue, :red)
    g = complete_digraph(5)
    f = Figure(size = (800, 800))
    graphplot(f[1,1], g, layout = Shell(), nlabels = nodes, edge_width = 3abs.(vals), edge_color = cols)

    return f
end

# visualize an aggregate from a vector of interaction matrices
function plot_densities(mat, nodes)
    mm = stack(mat)

    f = Figure(size = (900, 1200))
    for i in 1:5, j in 1:5
        if i != j
            col = sum(mm[i, j, :]) < 0 ? :red : :blue
            ax = Axis(f[i,j]; 
                title="$(nodes[i]) on $(nodes[j])", 
                limits=(-1, 1, nothing, nothing)
            )
            GLMakie.hist!(ax, mm[i,j,:], 
                bins=200, color=(col, 0.3), strokearound=true, strokewidth=1, strokecolor=col
            )
        end
    end
    f
end

# Metadata

function interaction_matrix(bodymass;
    competition_to_predation_max = 0.1
)
    # Create interaction matrix
    cats = 1
    black_rat = 2
    norway_rat = 3
    mouse = 4
    pig = 5
    rodents = [black_rat, norway_rat, mouse]

    interactions = zeros(5, 5)

    # cat/rodent
    # relative rankings of the effect of cats on an animal, low to high
    cat_order = [norway_rat, black_rat, mouse]
    interactions[cats, rodents] = -sort(rand(3))[sortperm(cat_order)]

    interactions[rodents, cats] .= (-interactions[cats, i] * sqrt(bodymass[i]) / sqrt(bodymass[cats]) for i in rodents)

    # cat/pig
    # Small negative interactions
    smallest_interaction = minimum(abs.(interactions[cats, black_rat:mouse]))
    # Pigs and cats are randomly less than half of the
    # smallest interactions of the other species
    interactions[cats, pig] = -rand() * 0.5smallest_interaction
    interactions[pig, cats] = -rand() * 0.5smallest_interaction

    # rodents
    ninteractions = 6
    low = 0
    high = competition_to_predation_max * abs(interactions[cats, black_rat])
    max_multiplier = 2
    randvals = -rejection_range(ninteractions, low, high, max_multiplier)
    for i in black_rat:mouse, j in black_rat:mouse
        if i != j
            interactions[i, j] = pop!(randvals) * sqrt(bodymass[i]) / sqrt(bodymass[j])
        end
    end

    # pig/rodent
    ninteractions = 6
    low = 0
    high = minimum(abs(x) for x in interactions if x!= 0)
    max_multiplier = 2
    randpigs = -rejection_range(ninteractions, low, high, max_multiplier)
    for i in black_rat:mouse
        interactions[i, 5] = pop!(randpigs)
        interactions[5, i] = pop!(randpigs)
    end

    return interactions
end

pred_names = ["cat", "black_rat", "norway_rat", "mouse", "pig"]
# run(`libreoffice tables/animals.csv`)
pred_df = filter(CSV.read("tables/animals.csv", DataFrame)) do row
    row.name in pred_names
end

bodymass = map(identity, pred_df.mass)

# one matrix
ints = interaction_matrix(bodymass)
plot_interactions(ints, pred_names)

# pred_pops = NamedVector{Tuple(map(Symbol, pred_names)),5}(rand(5)) .* 10
# pred_interactions = SMatrix{5,5}(interaction_matrix(bodymass))
# pred_interactions * pred_pops

# multiple
mat = map(_-> interaction_matrix(bodymass), 1:1000)
plot_densities(mat, pred_names)
map(1:100) do _
    interaction_matrix(bodymass)
end

using Unitful
using Optimization
using OptimizationOptimJL

function rodent_func(x, p)
    (; k, rt, nsteps, years) = p
    yield = x[1]
    N = k
    total_taken = zero(N)
    for i in 1:years
        for i in 1:nsteps
            N = (N .* k) ./ (N .+ (k .- N) .* exp.(.-(rt / (nsteps * u"yr^-1"))))
            taken = N * -yield
            N += taken 
            total_taken += taken
        end
    end
    return (taken=total_taken / years, final_pop=N)
end

rodent_take_unitless(x, p) = ustrip(rodent_take(x, p))
rodent_take(x, p) = rodent_func(x, p)[1]
rodent_pop(x, p) = rodent_func(x, p)[2]

rodent_df = pred_df[2:4, :]
rodent_names = map(Symbol, rodent_df.name)
R = NamedVector{Tuple(rodent_names),3}

cat_energy_intake = 2131u"kJ/d"
assimilation_efficiency = 0.84
rat_energy_content = 6.24u"kJ/g" 
rodent_mass = R(rodent_df.mass) .* u"g"
assimilated_energy_per_individual = rodent_mass .* rat_energy_content .* assimilation_efficiency
individuals_per_cat = cat_energy_intake ./ assimilated_energy_per_individual

rt = R(rodent_df.rmax) .* u"yr^-1"
k = R(rodent_df.carrycap) .* u"ha^-1"

# nsteps = 365
# 12 steps (months) is a close enough approximation of continuous
nsteps = 12
years = 10
fraction_eaten = 0.72 # McGregor 2015

# map([12, 365]) do nsteps
rodent_params = map(R(rodent_names)) do rodent
    (; k=k[rodent], rt=rt[rodent], nsteps, years)
end
p = first(rodent_params)
optimal_takes = map(rodent_params) do p
    x0 = [0.1]
    prob = OptimizationProblem(rodent_take_unitless, x0, p)
    rodent_take(x0, p)
    sol = solve(prob, NelderMead())
    sol.u
end

optimal_pops = map(rodent_pop, optimal_takes, rodent_params)
optimal_caught = .-map(rodent_take, optimal_takes, rodent_params) ./ u"yr"
max_supported_cats = uconvert.(u"km^-2", optimal_caught ./ individuals_per_cat .* fraction_eaten)
# max_supported_cats = uconvert.(u"ha^-1", optimal_caught ./ individuals_per_cat)

killer_pops = map(rodent_func, optimal_takes .* 2, rodent_params) |> NamedTuple |> pairs
killer_caught = .-map(rodent_take, optimal_takes .* 2, rodent_params) ./ u"yr"
killer_supported_cats = uconvert.(u"km^-2", killer_caught ./ individuals_per_cat .* fraction_eaten)

optimal_takes ./ 3 .* nsteps
realistic_pops = map(rodent_pop, optimal_takes ./ 3, rodent_params)
realistic_caught = .-map(rodent_take, optimal_takes ./ 3, rodent_params) ./ u"yr"
realistic_supported_cats = uconvert.(u"km^-2", realistic_caught ./ individuals_per_cat .* fraction_eaten)
realistic_supported_cats = uconvert.(u"ha^-1", realistic_caught ./ individuals_per_cat .* fraction_eaten)
# end

using Distributions
cat_preference = LogNormal(log(50), 1); 
p = Makie.plot(cat_preference; axis=(;yzoomlock=true))
xlims!(p.axis, (0, 300))
for x in 1:0.2:4.3
    Makie.plot!(p.axis, LogNormal(log(40), x^2/10))
    sleep(0.1)
end

cat_mean_prey_sizes = [11.0, 33.0, 17.9, 4.3, 15.5, 23.6, 24.0, 31.2, 26.5, 23.3, 9.3, 26.2, 29.0, 21.2, 349.7, 3.64, 21.8, 35.9, 13.6, 8.2, 228.6, 249.0, 42.7, 16.5, 15.5, 60.7, 32.8, 184.9, 34.1, 36.8, 7.4, 102.0, 241.2, 38.5, 13.3, 26.4, 7.6, 24.3, 123.2, 39.2, 1.77, 76.7, 13.5, 18.3, 30.3, 34.4]
sort!(cat_mean_prey_sizes)
mean(cat_mean_prey_sizes)
std(cat_mean_prey_sizes)
sort!(rand(cat_preference, length(cat_mean_prey_sizes)))

# Wilson et al 2007
black_rat_mass = 20:40:220
black_rat_frequency = [0.08, 0.20, 0.26, 0.31, 0.14, 0.01]
mean_black_rat_mass = sum(black_rat_frequency .* black_rat_mass)
Makie.plot(LogNormal(log(mean_black_rat_mass)))
vals = map(1:500) do x
    pdf(SkewNormal(mean_black_rat_mass, 30, -10), x)
end
Makie.lines(vals)

x0 = [20, -10]
p = (; cat_mean_prey_sizes)
prey_size(x0, p)
prob = OptimizationProblem(prey_size, x0, p)
sol = solve(prob, NelderMead(); )

function prey_size(x, p)
    cat_preference = LogNormal(log(50), x[1]); 
    sum(1:10000) do _
        sum(abs.(sort!(rand(cat_preference, length(p.cat_mean_prey_sizes))) .- p.cat_mean_prey_sizes))
    end
end

x0 = [0.1]
p = (; cat_mean_prey_sizes)
prey_size(x0, p)
prob = OptimizationProblem(prey_size, x0, p)
sol = solve(prob, NelderMead(); )
cat_preference = LogNormal(log(50), 0.79); 
Makie.hist!(cat_mean_prey_sizes)
p = Makie.plot(cat_preference)
xlims!(p.axis, (0, 300))

rand(cat_preference)
cdf(cat_preference, 200) - cdf(cat_preference, 400)
cdf(cat_preference, 50) - cdf(cat_preference, 20) 


# Write a model of rodent/rabbit preference given populations
function catch((; black_rat, norway_rat, mouse))

end

#= 
Another paper using this simulation:

Mesopredator release of rats can only happen when there is a high-density
primary food-source for cats *other than rodents that compete with each other*. 
This is likely limited to urban feeding/refuse, 
high density sea bird (and similar) colonies, or rabbits.

Demonstration: island with or without rabbits and only rats
- grasslands with rabbits will show meopredator release
- forests will not show mesopredator release
=#
