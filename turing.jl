# second example, with a vector of parameters
using Turing, Statistics, Distributions, Random, LinearAlgebra, StatsPlots
using AdvancedMH
using Plots
using LinearAlgebra

cellsize = 1
elev_center = 1
slopefactor = -3.5
distfactor = 0.05
function _slope_speed(e, d=10)
    slope = (e - elev_center) / (cellsize * d)
    slope = 
end

Plots.plot(Base.Fix2(_slope_speed, 2), 0:0.1:10)

logit(x) = one(x)/(one(x) + exp(-x))

function predict_extinctions(ruleset, init; tspan, kw...)
    output = TransformedOutput(init; tspan, kw...) do data
        map(>(0), sum(data.species))
    end
    sim!(output, ruleset)
    # Return NamedVector of years to extinction
    return sum(sim) .+ first(tspan) .- 1
end

function define_threats(human, cat, rat, habitat)
    # define a threat rule here
    # probably should not be a new anonymous function
    # or we will get recompilation problems
end

@model function extinction_model(sim, traits, priors, extinct_dates)
    ##################################33
    # Priors
    
    # Predator/Human level parameters for each trait
    cat_target_mass ~ Normal(priors.mass.target.cat...)
    rat_target_mass ~ Normal(priors.mass.target.rat...)
    human_target_mass ~ Normal(priors.mass.target.human...)
    target_masses = (cat_target_mass, rat_target_mass, human_target_mass)

    cat_target_slope ~ Beta(priors.mass.target.cat...)
    rat_target_slope ~ Beta(priors.mass.target.rat...)
    human_target_slope ~ Beta(priors.mass.target.human...)
    mass_slopes = (cat_target_slope, rat_target_slope, human_target_slope) 

    cat_mass_intensity ~ LogNormal(priors.mass.target.cat...)
    rat_mass_intensity ~ LogNormal(priors.mass.target.rat...)
    human_mass_intensity ~ LogNormal(priors.mass.target.human...)
    max_mass_intensity = (cat_target_slope, rat_target_slope, human_target_slope) 

    pred_groundnesting_effects ~ filldist(LogNormal(0, 2), 3)
    pred_flight_effects ~ filldist(LogNormal(0, 2), 3)
    # pred_class_effects .~ Normal.(priors.class, 3)
    # Human edibility effect
    edibility_effect ~ Normal(priors.edibility)

    # How much native habitat is favoured
    habitat .~ Beta.(priors.habitat, 1.0)

    ##################################33
    # Combine binary parameters with values
    
    spec_flight_effects = map(.*, traits.flight, pred_flight_effects)
    spec_groundnesting_effects = map(.*, traits.groundnesting, pred_groundnesting_effects)

    # Species-level parameters
    spec_combined_effects = map(traits.log_mass, spec_flight_effects, spec_groundnesting_effects) do mass, flight, goundnest
        map(target_mass, mass_slopes, max_mass_intensity) do target, slop, intensity
            distance = abs(mass - target)
            slope^distance * intensity * flight * groundnest
        end
    end

    spec_edibility_effects = map(.*, traits.edibility, edibility_effect)
    human = spec_combined_effects.human .+ edibility_effect
    # Define the threat rule with these values
    threatrule = define_threats(human, spec_combined_effects, habitat)
    # build a new ruleset
    ruleset = DynamicGrids.StaticRuleset(threatrule, otherrules...)
    # Run the simulation `nreplicates` times
    mean_predictions = mean(1:nreplicates) do
        predict_extinctions(sim, encounter_probabilities)
    end
    return extinct_dates .~ Normal.(mean_predictions, 10)
end

Plots.plot(Beta(4, 6))

mus_native = deepcopy(obs.mus.native)
mus_native_filled = mapcols(deepcopy(mus_native)) do col
    if eltype(col) == Bool
        l = findlast(col)
        isnothing(l) && return col
        range = 1:l
        col[range] .= true
    end
    return col
end
species_names = Symbol.(names(mus_native_filled)[2:20])
obs = cat(getindex.(Ref(mus_native), :, species_names)...; dims=Dim{:species}(species_names))
preds = cat(getindex.(Ref(mus_native_filled), :, species_names)...; dims=Y(species_names))
humans = human_pop_timeline.mus[Near(lookup(mus_native.Raphus_cucullatus, Ti))]

mascarene_species_csv = "/home/raf/PhD/Mascarenes/MauritiusExtinctions/mascarine_species.csv"
mascarene_species = CSV.read(mascarine_species_csv, DataFrame)

using StatsPlots, Distributions
using Distributions: Normal
Plots.plot(Normal(log(100), 10))

priors = ( 
    mass = (
        target_mass = (
           rat = (mean=log(10), var=1),
           cat = (mean=log(41.19), var=1), # (PEARRE & MAASS, 1998 - mean prey size at ~20° latitude)
           human = (mean=log(100), var=1), 
        ),
        effect = (
           rat = (mean=10, var=10),
           cat = (mean=30, var=10), # (PEARRE & MAASS, 1998 - mean prey size at ~20° latitude)
           human = (mean=30, var=10), 
        ),
    )
)

exp(1)
vals = rand(LogNormal(10, 1), 100000)
Plots.histogram(log.(vals))
        
using BenchmarkTools

using FillArrays
using Turing, Distributions
using MCMCChains, Plots, StatsPlots
using StatsFuns: logistic
using MLDataUtils: shuffleobs, stratifiedobs, rescale!
using Random

# @views function main()
θ = [2.0; 3.0]
n = 1000 # sample size
S = 1000 # number of simulation draws
# the data generating process
function dgp(θ, n)
    [rand(Exponential(θ[1]), n) rand(Poisson(θ[2]), n)]
end

# summary statistics for estimation
function moments(y)
    sqrt(n) .* [mean(y, dims=1)[:]; std(y, dims=1)[:]]
end

y = dgp(θ, n)
z = moments(y)

@model function abc(z)
    zs = zeros(S, size(z,1))
    # create the prior: the product of the following array of marginal priors
    θ ~ arraydist([LogNormal(1.,1.); LogNormal(1.,1.)])

    # sample from the model, at the trial parameter value, and compute statistics
    @inbounds for i = 1:S
        y .= dgp(θ, n)
        zs[i, :] .= moments(y) # simulated summary statistics
    end
    # the asymptotic Gaussian distribution of the statistics
    m = mean(zs, dims=1)[:]
    Σ = Symmetric((1.0 + 1/S) * cov(zs))
    z ~ MvNormal(m, Σ)
end

# sample chains, 4 chains of length 1000, with 100 burnins dropped
len = 1000
burnin = 100
chain = sample(abc(z), MH(:θ => AdvancedMH.RandomWalkProposal(MvNormal(zeros(2), 0.25 * I))), MCMCThreads(), len+burnin, 4)

using AdvancedMH
chain = chain[burnin+1:end, :, :]
display(chain)
p = Plots.plot(chain)
display(p)

?AdvancedMH.RandomWalkProposal

using Distributions, Turing
y = rand(filldist(Normal(2.5, 0.5), 100))
using StatsPlots

Plots.plot(truncated(Normal(0, 1); lower = 0.0))
function mymodel(x)
    sum(rand(x))
end

using PDMats, LinearAlgebra
@btime f(zs)
using DimensionalData
zs = rand((S), (2000))
p = PDMat(Cholesky(LowerTriangular(cov(zs; dims=1))))
(1 + 1/S) * cov(zs)
function f(zs)
    S = 10
    zs = rand(S, 200)
    Σ = PDMat(Cholesky(LowerTriangular((1.0 + 1/S) * cov(zs))))
    m = vec(mean(zs, dims=1))
    MvNormal(m, Σ)
end

@macroexpand
@model function test_model(y)
    λ ~ Poisson(4)
    error ~ truncated(Normal(0, 1); lower = 0.0)
    y_hat = [mymodel(λ) for i in 1:100]
    y ~ MvNormal(y_hat, error)
end
logpdf(distr, y)

model = test_model(y)
sample(model, PG(10), 100)

vals = [mymodel(2) for i in 1:100]

MvNormal(y_hat, error)

using Turing, Distributions
using StatsPlots

@model function gdemo(y)
    s² ~ InverseGamma(1, 1)
    m ~ Normal(0.5, sqrt(s²))
    for i in eachindex(y)
        y[i] ~ Normal(m, sqrt(s²))
    end
end
a = rand(1000)
chn = sample(gdemo(a), NUTS(), 1000)
var(a)
mean(a)
plot(chn)
s² = InverseGamma(2, 3)
m = Normal(0.5, 1)
x_bar = mean(a)
N = length(a)

mean_exp = (m.σ * m.μ + N * x_bar) / (m.σ + N)

plot(InverseGamma(1, 1); xlims=(0, 3))
plot(Normal(1, 20))
Normal(m, sqrt(s²))
return y ~ Normal(m, sqrt(s²))

# Import libraries.
using Turing, StatsPlots, Random

# Set the true probability of heads in a coin.
p_true = 0.8

# Iterate from having seen 0 observations to 100 observations.
Ns = 0:100

# Draw data from a Bernoulli distribution, i.e. draw heads or tails.
Random.seed!(12)
data = rand(Bernoulli(p_true), last(Ns))

using Turing, Statistics, Distributions, Random, LinearAlgebra, StatsPlots
using AdvancedMH
# Declare our Turing model.
@model function coinflip(y)
    # Our prior belief about the probability of heads in a coin.
    p ~ Beta(1, 1)

    # The number of observations.
    N = length(y)
    # Heads or tails of a coin are drawn from a Bernoulli distribution.
    y .~ Bernoulli.(p)
end

# Settings of the Hamiltonian Monte Carlo (HMC) sampler.
iterations = 1000
ϵ = 0.05
τ = 10

# Start sampling.
chain = sample(coinflip(data), HMC(ϵ, τ), iterations)
len = 1000
burnin = 100
chain = sample(coinflip(data), 
    MH(:p => AdvancedMH.RandomWalkProposal(MvNormal(zeros(2), 0.25*I))),
    MCMCThreads(), len+burnin, 4
)
chain = chain[burnin+1:end,:,:]
display(chain)

# Plot a summary of the sampling process for the parameter p, i.e. the probability of heads in a coin.
histogram(chain[:p])


@model function gdemo(y)
    s ~ InverseGamma(2,3)
    m ~ Normal(1.7, sqrt(s))
    y .~ Normal.(m, sqrt(s))
end

# Use a static proposal for s and random walk with proposal
# standard deviation of 0.25 for m.
data = rand(1000) .+ 1.2
mean(data)
chain = sample(
  gdemo(data),
  MH(
      :s => AdvancedMH.StaticProposal(InverseGamma(2,3)),
      :m => AdvancedMH.RandomWalkProposal(Normal(0, 0.25))
  ),
  10_000
)
mean(chain)


# Import Turing.
using Turing
# Package for loading the data set.
using RDatasets
# Package for visualization.
using StatsPlots
# Functionality for splitting the data.
using MLUtils: splitobs
# Functionality for constructing arrays with identical elements efficiently.
using FillArrays
# Functionality for normalizing the data and evaluating the model predictions.
using StatsBase
# Functionality for working with scaled identity matrices.
using LinearAlgebra
# Set a seed for reproducibility.
using Random
Random.seed!(0);
using Distributions: Normal

# Load the dataset.
data = RDatasets.dataset("datasets", "mtcars")

# Show the first six rows of the dataset.
first(data, 6)

# Remove the model column.
select!(data, Not(:Model))

# Split our dataset 70%/30% into training/test sets.
trainset, testset = map(DataFrame, splitobs(data; at=0.7, shuffle=true))

# Turing requires data in matrix form.
target = :MPG
train = Matrix(select(trainset, Not(target)))
test = Matrix(select(testset, Not(target)))
train_target = trainset[:, target]
test_target = testset[:, target]

# Standardize the features.
dt_features = fit(ZScoreTransform, train; dims=1)
StatsBase.transform!(dt_features, train)
StatsBase.transform!(dt_features, test)

# Standardize the targets.
dt_targets = fit(ZScoreTransform, train_target)
StatsBase.transform!(dt_targets, train_target)
StatsBase.transform!(dt_targets, test_target);

# Bayesian linear regression.
@macroexpand @model function linear_regression(x, y)
    # Set variance prior.
    σ² ~ truncated(Normal(0, 100); lower=0)

    # Set intercept prior.
    intercept ~ Normal(0, sqrt(3))

    # Set the priors on our coefficients.
    nfeatures = size(x, 2)
    coefficients ~ MvNormal(Zeros(nfeatures), 10.0 * I)

    # Calculate all the mu terms.
    mu = intercept .+ x * coefficients
    return y ~ MvNormal(mu, σ² * I)
end
intercept = rand(Normal(0, sqrt(3)))
coef = rand(MvNormal(Zeros(10), 10.0 * I))
intercept .+ train * coef
[1 2; 3 4; 5 6] * [10, 20]
train_target
Plots.plot(truncated(Normal(0, 100); lower=0))

model = linear_regression(train, train_target)
chain = sample(model, NUTS(), 500)
u = I * 10

#Import Turing, Distributions and DataFrames
using Turing, Distributions, DataFrames, Distributed

# Import MCMCChain, Plots, and StatsPlots for visualizations and diagnostics.
using MCMCChains, Plots, StatsPlots

# Set a seed for reproducibility.
using Random
Random.seed!(12);

theta_noalcohol_meds = 1    # no alcohol, took medicine
theta_noalcohol_meds = 1    # no alcohol, took medicine
theta_alcohol_meds = 3      # alcohol, took medicine
theta_noalcohol_nomeds = 6  # no alcohol, no medicine
theta_alcohol_nomeds = 36   # alcohol, no medicine

# no of samples for each of the above cases
q = 100

#Generate data from different Poisson distributions
noalcohol_meds = Poisson(theta_noalcohol_meds)
alcohol_meds = Poisson(theta_alcohol_meds)
noalcohol_nomeds = Poisson(theta_noalcohol_nomeds)
alcohol_nomeds = Poisson(theta_alcohol_nomeds)

nsneeze_data = vcat(
    rand(noalcohol_meds, q),
    rand(alcohol_meds, q),
    rand(noalcohol_nomeds, q),
    rand(alcohol_nomeds, q),
)
alcohol_data = vcat(zeros(q), ones(q), zeros(q), ones(q))
meds_data = vcat(zeros(q), zeros(q), ones(q), ones(q))

df = DataFrame(;
    nsneeze=nsneeze_data,
    alcohol_taken=alcohol_data,
    nomeds_taken=meds_data,
    product_alcohol_meds=meds_data .* alcohol_data,
)
df[sample(1:nrow(df), 5; replace=false), :]

#Data Plotting

p1 = Plots.histogram(
    df[(df[:, :alcohol_taken] .== 0) .& (df[:, :nomeds_taken] .== 0), 1];
    title="no_alcohol+meds",
)
p2 = Plots.histogram(
    (df[(df[:, :alcohol_taken] .== 1) .& (df[:, :nomeds_taken] .== 0), 1]);
    title="alcohol+meds",
)
p3 = Plots.histogram(
    (df[(df[:, :alcohol_taken] .== 0) .& (df[:, :nomeds_taken] .== 1), 1]);
    title="no_alcohol+no_meds",
)
p4 = Plots.histogram(
    (df[(df[:, :alcohol_taken] .== 1) .& (df[:, :nomeds_taken] .== 1), 1]);
    title="alcohol+no_meds",
)
Plots.plot(p1, p2, p3, p4; layout=(2, 2), legend=false)

# Convert the DataFrame object to matrices.
df
data = Bool.(Matrix(df[:, [:alcohol_taken, :nomeds_taken, :product_alcohol_meds]]))
data_labels = df[:, :nsneeze]
data

# Bayesian poisson regression (LR)
@model function poisson_regression(x, y, n, σ²)
    b0 ~ Normal(0, σ²)
    b1 ~ Normal(0, σ²)
    b2 ~ Normal(0, σ²)
    b3 ~ Normal(0, σ²)
    for i in 1:n
        theta = b0 + b1 * x[i, 1] + b2 * x[i, 2] + b3 * x[i, 3]
        y[i] ~ Poisson(exp(theta))
    end
end;

# Retrieve the number of observations.
n, _ = size(data)

# Sample using NUTS.

num_chains = 4
m = poisson_regression(data, data_labels, n, 10)
chain = sample(m, NUTS(), MCMCThreads(), 2000, num_chains; discard_adapt=false)
chain


# Taking the first chain
c1 = chain[:, :, 1]

# Calculating the exponentiated means
b0_exp = exp(mean(c1[:b0]))
b1_exp = exp(mean(c1[:b1]))
b2_exp = exp(mean(c1[:b2]))
b3_exp = exp(mean(c1[:b3]))

print("The exponent of the meaned values of the weights (or coefficients are): \n")
println("b0: ", b0_exp)
println("b1: ", b1_exp)
println("b2: ", b2_exp)
println("b3: ", b3_exp)
print("The posterior distributions obtained after sampling can be visualised as :\n")

plot(chain)
