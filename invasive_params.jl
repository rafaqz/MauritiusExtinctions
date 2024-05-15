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
    norway_rat = 2
    black_rat = 3
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
using Distributions
using Optimization
using OptimizationOptimJL


function rodent_func(x, p)
    replicates = map(1:500) do i
        (; k, rt, nsteps, years) = p
        N = k
        total_taken = zero(N)
        for i in 1:years
            for i in 1:nsteps
                yield = rand(LogNormal(log(x[1]), p.std))
                N = (N .* k) ./ (N .+ (k .- N) .* exp.(.-(rt / (nsteps * u"yr^-1"))))
                # Cant take the whole population
                taken = max(-N + oneunit(N) * 0.1, N * -yield)
                @assert taken <= zero(taken) "taken: $taken N: $N yield: $yield"
                @assert taken >= -N "taken: $taken N: $N yield: $yield"
                N += taken
                total_taken += taken
            end
        end
        return (taken=total_taken / years, final_pop=N)
    end
    return (taken=mean(first, replicates), final_pop=mean(last, replicates))
end

pred_names = ["cat", "black_rat", "norway_rat", "mouse", "pig"]
# run(`libreoffice tables/animals.csv`)
pred_df = filter(CSV.read("tables/animals.csv", DataFrame)) do row
    row.name in pred_names
end

rodent_take_unitless(x, p) = ustrip(rodent_take(x, p))
rodent_take(x, p) = rodent_func(x, p)[1]
rodent_pop(x, p) = rodent_func(x, p)[2]

rodent_df = pred_df[2:4, :]
rodent_names = map(Symbol, rodent_df.name)
rodent_labels = titlecase.(replace.(rodent_df.name, ("_" => " ",)))
R = NamedVector{Tuple(rodent_names),3}

cat_energy_intake = 2131u"kJ/d"
assimilation_efficiency = 0.84
rat_energy_content = 6.24u"kJ/g" 
# rodent_mass = R(rodent_df.mass) .* u"g"
# Use mean from actual sizes taken for Norway rats, and guess the others
# For mouse we just use the mean size.
rodent_mass = NamedVector(norway_rat=110, black_rat=90, mouse=27) .* u"g"
assimilated_energy_per_individual = rodent_mass .* rat_energy_content .* assimilation_efficiency
individuals_per_cat = cat_energy_intake ./ assimilated_energy_per_individual

rt = R(rodent_df.rmax) .* u"yr^-1"
k = R(rodent_df.carrycap) .* u"ha^-1"

# nsteps = 365
# 12 steps (months) is a close enough approximation of continuous
nsteps = 12
years = 10
fraction_eaten = 0.72 # McGregor 2015

sds = (0.0:0.1:1.5).^2
stochastic_rates = map(sds) do std
    rodent_params = map(R(rodent_names)) do rodent
        (; k=k[rodent], rt=rt[rodent], nsteps, years, std)
    end
    p = first(rodent_params)
    optimal_takes = map(rodent_params) do p
        x0 = [0.1]
        prob = OptimizationProblem(rodent_take_unitless, x0, p; lb=[0.0], ub=[1.0])
        rodent_take(x0, p)
        sol = solve(prob, SAMIN(); maxiters=10000)
        sol.u
    end
    optimal_pops = map(rodent_pop, optimal_takes, rodent_params)
    optimal_caught = map(rodent_take, optimal_takes, rodent_params) ./ u"yr" .* -1
    max_supported_cats = uconvert.(u"km^-2", optimal_caught ./ individuals_per_cat .* fraction_eaten)
    @show optimal_pops optimal_caught max_supported_cats
    takes = NamedTuple{map(k -> Symbol(:taken_, k), propertynames(optimal_takes))}(values(optimal_takes))
    cats = NamedTuple{map(k -> Symbol(:cat_, k), propertynames(max_supported_cats))}(values(max_supported_cats))
    (; std, map(first, takes)..., cats..., NamedTuple(optimal_pops)...)
end |> DataFrame

colors = [:red, :lightblue, :yellow]
x = stochastic_rates.std
fig = Figure(; title="Effects of monthly stochasticity on optimal hunting yields", size=(600, 800));
labels = ["Mean fraction hunted per month", "Supported cats per km^2", "Rodents per hectare"]
axs = map(1:3, labels) do i, ylabel
    ax = if i == 3
        Axis(fig[i, 1]; 
            xlabel="Scale of log-normal hunting stochasticity", ylabel,
        )
    else
        Axis(fig[i, 1]; ylabel, xticksvisible=false, xticklabelsvisible=false)
    end
    # hidespines!(ax)
    ax
end
rodent_plots = map(1:3) do i
    color = colors[i]
    label = rodent_labels[i]
    p1 = Makie.lines!(axs[1], x, ustrip.(stochastic_rates[:, i+1]); 
        label, color
    )
    p2 = Makie.lines!(axs[2], x, ustrip.(stochastic_rates[:, i+4]); 
        label, color
    )
    p3 = Makie.lines!(axs[3], x, ustrip.(stochastic_rates[:, i+7]); 
        label, color
    )
    (p1, p2, p3)
end
axislegend(axs[1])
save("images/cat_hunting_stochasticty.png", fig)

using Distributions
using NLopt

# Childs 1986
max_rat_size = 100:100:700
center_rat_size = 50:100:650
max_rat_25_size = 25:25:200
center_rat_25_size = 12.5:25:187.5
killed_rats_25_childs = [0, 8, 8, 3, 1, 1, 0, 1]
killed_rats_childs = [19, 4, 0, 0, 0, 0, 0]
trapped_rats_childs = [5, 5, 9, 29, 44, 14, 3]
sum(center_rat_25_size .* killed_rats_25_childs) ./ sum(killed_rats_25_childs)

# Glass 2009
trapped_rats_glass = [18, 38, 49, 167, 175, 98, 0]
killed_rats_glass = [15, 13, 4, 2, 1, 0, 0]
center_rat_size .* trapped_rats_glass
center_rat_size .* killed_rats_glass
cat_preference_glass = LogNormal(log(60), 1.1); 
rat_pdfs_glass = pdf.((cat_preference_glass,), center_rat_size) .* trapped_rats_glass
killed_rats_glass ./ sum(killed_rats_glass) 
rat_pdfs_glass ./ sum(rat_pdfs_glass)

mean_glass = sum(center_rat_size .* killed_rats_glass) ./ sum(killed_rats_glass)
mean_childs = sum(center_rat_25_size .* killed_rats_25_childs) ./ sum(killed_rats_25_childs)
n_glass = sum(killed_rats_glass) 
n_childs = sum(killed_rats_25_childs)
n_total = n_glass + n_childs 
frac_glass = n_glass / n_total
frac_childs = n_childs / n_total
mean_caught_norway_rat_mass = mean_glass .* frac_glass + mean_childs .* frac_childs

trapped_rats = trapped_rats_glass .+ trapped_rats_childs
killed_rats = killed_rats_glass .+ killed_rats_childs

function get_pref(x, p)
    μ, σ = x
    cat_preference = LogNormal(log(μ), σ); 
    rat_pdfs = pdf.((cat_preference,), p.center_rat_size) .* p.trapped_rats
    killed_means = p.killed_rats ./ sum(p.killed_rats)
    predicted_means = rat_pdfs ./ sum(rat_pdfs)
    sum((killed_means .- predicted_means) .^ 2)
end

x = [90.0, 0.6]
p = (; trapped_rats, killed_rats, center_rat_size)
prob = OptimizationProblem(get_pref, x, p; lb=[0.0, 0.0], ub=[500.0, 10.0])
get_pref(x, p)
sol = solve(prob, NLopt.LN_NELDERMEAD())
x = sol.u
get_pref(x, p)
μ, σ = x
cat_preference = LogNormal(log(μ), σ)
Makie.plot(cat_preference)


# From the diet review paper
# these are means prey size accross whole studies
using Statistics
cat_mean_prey_sizes = [11.0, 33.0, 17.9, 4.3, 15.5, 23.6, 24.0, 31.2, 26.5, 23.3, 9.3, 26.2, 29.0, 21.2, 349.7, 3.64, 21.8, 35.9, 13.6, 8.2, 228.6, 249.0, 42.7, 16.5, 15.5, 60.7, 32.8, 184.9, 34.1, 36.8, 7.4, 102.0, 241.2, 38.5, 13.3, 26.4, 7.6, 24.3, 123.2, 39.2, 1.77, 76.7, 13.5, 18.3, 30.3, 34.4]
sort!(cat_mean_prey_sizes)
mean(cat_mean_prey_sizes)
Statistics.std(cat_mean_prey_sizes)
sort!(rand(cat_preference, length(cat_mean_prey_sizes)))

x0 = [20.0, -10.0]
p = (; cat_mean_prey_sizes)
prey_size(x0, p)
prob = OptimizationProblem(prey_size, x0, p)
sol = solve(prob, NelderMead())

x0 = [0.1]
p = (; cat_mean_prey_sizes)
prey_size(x0, p)
prob = OptimizationProblem(prey_size, x0, p)
sol = solve(prob, NelderMead(); )
cat_preference = LogNormal(log(50), 0.8); 
Makie.hist!(cat_mean_prey_sizes)
p = Makie.plot(cat_preference)
xlims!(p.axis, (0, 300))

# Wilson et al 2007
black_rat_max_mass = 40:40:240
black_rat_center_mass = 20:40:220
black_rat_trap_rate = [0.08, 0.20, 0.26, 0.31, 0.14, 0.01]
mean_black_rat_mass = sum(black_rat_trap_rate .* black_rat_center_mass)
black_rat_mass_distribution = black_rat_trap_rate .* black_rat_center_mass ./ mean_black_rat_mass

# Glass 1988
norway_rat_max_mass = 100:100:500
norway_rat_center_mass = norway_rat_max_mass .- 50
male_norway_rat_trap_rate = [0.07, 0.3, 0.43, 0.16, 0.04]
female_norway_rat_trap_rate = [0.14, 0.33, 0.34, 0.19, 0.0]
norway_rat_trap_rate = (male_norway_rat_trap_rate .+ female_norway_rat_trap_rate) ./ 2 
mean_norway_rat_mass = sum(norway_rat_trap_rate .* norway_rat_center_mass)
norway_rat_mass_distribution = norway_rat_trap_rate .* norway_rat_center_mass ./ mean_norway_rat_mass

# made up from shifted rat samples
mouse_max_mass = 10:10:60
mouse_center_mass = 5:10:55
mouse_trap_rate = [0.08, 0.20, 0.26, 0.31, 0.14, 0.01]
mean_mouse_mass = sum(mouse_trap_rate .* mouse_center_mass)
mouse_mass_distribution = (mouse_trap_rate .* mouse_center_mass) ./ mean_mouse_mass
mouse_mass_yields = mouse_rates .* mouse_mass_distribution

norway_rat_pdfs = pdf.((cat_preference,), norway_rat_center_mass)
norway_rat_trap_rate
norway_rat_catch_rates = norway_rat_pdfs .* norway_rat_trap_rate
norway_rat_normalised_rates = norway_rat_catch_rates ./ sum(norway_rat_rates)
norway_rat_normalised_mass_yield = sum(norway_rat_normalised_rates .* norway_rat_mass_distribution .* mean_norway_rat_mass)
norway_rat_mass_yields = sum(norway_rat_catch_rates .* norway_rat_mass_distribution)
norway_rat_predation_rate = sum(norway_rat_catch_rates)

black_rat_pdfs = pdf.((cat_preference,), black_rat_center_mass)
black_rat_trap_rate
black_rat_catch_rates = black_rat_pdfs .* black_rat_trap_rate
black_rat_normalised_rates = black_rat_catch_rates ./ sum(black_rat_catch_rates)
black_rat_normalised_mass_yield = sum(black_rat_normalised_rates .* black_rat_mass_distribution .* mean_black_rat_mass)
black_rat_mass_yield = sum(black_rat_catch_rates .* black_rat_mass_distribution)
black_rat_predation_rate = sum(black_rat_catch_rates)

mouse_pdfs = pdf.((cat_preference,), mouse_center_mass)
mouse_trap_rate
mouse_catch_rates = mouse_pdfs .* mouse_trap_rate
mouse_normalised_rates = mouse_catch_rates ./ sum(mouse_catch_rates)
mouse_normalised_mass_yield = sum(mouse_normalised_rates .* mouse_mass_distribution .* mean_mouse_mass)
mouse_mass_yield = sum(mouse_catch_rates .* mouse_mass_distribution)
mouse_predation_rate = sum(mouse_catch_rates)

norway_rat_predation_rate / black_rat_predation_rate 
norway_rat_predation_rate / mouse_predation_rate 
black_rat_predation_rate / mouse_predation_rate 
means = (mean_norway_rat_mass, mean_black_rat_mass, mean_mouse_mass)


# norway_rat : black_rat : mouse
# Roughly ~ 1:3:10
#
# rodents are ~80% of cat diet and they do not switch to natives much (Harper thesis etc)

fig = Figure(; size=(600, 800));
ax1 = Axis(fig[1, 1]; ylabel="Cat mass preference")
ax2 = Axis(fig[1, 1]; ylabel="Rodent mass fractions", yaxisposition=:right)
# ax3 = Axis(fig[2, 1])
hidexdecorations!(ax2)
hidespines!(ax2)
linkxaxes!(ax1, ax2, ax3)
alpha = 0.7
b1 = Makie.lines!(ax2, norway_rat_center_mass, norway_rat_trap_rate; 
    color=(colors[1], alpha), label=rodent_labels[1] * " trapped sizes",
)
b2 = Makie.lines!(ax2, black_rat_center_mass, black_rat_trap_rate; 
    color=(colors[2], alpha), label=rodent_labels[2] * " trapped sizes",
)
b3 = Makie.lines!(ax2, mouse_center_mass, mouse_trap_rate; 
    color=(colors[3], alpha), label=rodent_labels[3] * " trapped sizes",
)
l = Makie.plot!(ax1, cat_preference; color=:black, label="Log-normal preference model")
d = Makie.density!(ax2, cat_mean_prey_sizes; color=:grey, label="Literature mean prey sizes")
axislegend(ax1; position=:rc)
axislegend(ax2)
xlims!(ax1, (0, 500))
xlims!(ax2, (0, 500))
# xlims!(ax3, (0, 500))
save("images/cat_rodent_predation.png", fig)

Makie.lines!(ax3, center_rat_size, killed_rats_glass ./ sum(killed_rats_glass))
Makie.barplot!(ax3, center_rat_size, trapped_rats_childs ./ sum(trapped_rats_childs))
Makie.lines!(ax3, center_rat_25_size, killed_rats_25_childs ./ sum(killed_rats_25_childs) .* 4)
yaxisright!(ax2)
display(fig)


# Write a model of rodent/rabbit preference given populations
function catch()

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
