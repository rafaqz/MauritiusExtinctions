# second example, with a vector of parameters
using Turing, Statistics, Distributions, Random, LinearAlgebra, StatsPlots
using AdvancedMH
using Plots
using LinearAlgebra

# Encounter probability is the inverse of 
# the square of the cost distance?
function encounter_probability(cost_distance)
    1 / cost_distance ^ 2
end

function obs_probability(ruleset, init, encounter_probabilities, visibility; kw...)
    output = TransformedOutput(init; kw...) do data
        encounter_probs = view(encounter_probabilities, Ti(Contains(timestep(data))))
        sum(splat(*), zip(data.species, encounter_probs, visibility))
    end
    return sim!(output, ruleset)
end

function define_threats(human, cat, rat, habitat)
    # define a threat rule here
    # probably should not be a new anonymous function
    # or we will get recompilation problems
end

@model function species_observations(ruleset, init, encounter_probabilities, visibility, species; kw...) where Keys
    # Observation visibility prior for each species
    visibility .~ Normal.(species.visibility, 1.0)
    # Need some distributions between zero and one...
    cat .~ Binomial.(species.susceptibility.cat, 1.0)
    rat .~ Binomial.(species.susceptibility.rat, 1.0)
    human .~ Binomial.(species.susceptibility.human, 1.0)
    # How much native habitat is favoured
    habitat .~ Binomial.(species.susceptibility.habitat, 1.0)
    # Define the threat rule with these values
    threatrule = define_threats(human, cat, rat, habitat)
    # build a new ruleset
    ruleset = DynamicGrids.StaticRuleset(threatrule, otherrules...)
    # Run the simulation
    for i = 1:S
        # obs is a DimArray over time of NamedVector for each species
        obs_predicted = obs_probability(ruleset, init, encounter_probabilities, visibility; kw...)
        sum(zip(obs_recorded, obs_predicted)) do (r::Bool, p::AbstractFloat)
            # some loss function 
        end
    end
end

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

@macroexpand @model function abc(z)
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

@macroexpand 
@model function test_model(y)
    λ ~ Poisson(4)
    error ~ truncated(Normal(0, 1); lower = 0.0)    
    y_hat = [mymodel(λ) for i in 1:100]    
    y ~ MvNormal(y_hat, error)
end

model = test_model(y)
sample(model, PG(10), 100)

vals = [mymodel(2) for i in 1:100]    

MvNormal(y_hat, error)


function test_model(__model__::DynamicPPL.Model, __varinfo__::DynamicPPL.AbstractVarInfo, __context__::AbstractPPL.AbstractContext, y; )
    @show y
    #= REPL[32]:1 =#
    begin
        #= REPL[32]:1 =#
        #= REPL[32]:2 =#
        begin
            var"##dist#381" = Poisson(4)
            var"##vn#378" = (DynamicPPL.resolve_varnames)((AbstractPPL.VarName){:λ}(), var"##dist#381")
            var"##isassumption#379" = begin
                    if (DynamicPPL.contextual_isassumption)(__context__, var"##vn#378")
                        if !((DynamicPPL.inargnames)(var"##vn#378", __model__)) || (DynamicPPL.inmissings)(var"##vn#378", __model__)
                            true
                        else
                            λ === missing
                        end
                    else
                        false
                    end
                end
            if (DynamicPPL.contextual_isfixed)(__context__, var"##vn#378")
                λ = (DynamicPPL.getfixed_nested)(__context__, var"##vn#378")
            elseif var"##isassumption#379"
                begin
                    (var"##value#382", __varinfo__) = (DynamicPPL.tilde_assume!!)(__context__, (DynamicPPL.unwrap_right_vn)((DynamicPPL.check_tilde_rhs)(var"##dist#381"), var"##vn#378")..., __varinfo__)
                    λ = var"##value#382"
                    var"##value#382"
                end
            else
                if !((DynamicPPL.inargnames)(var"##vn#378", __model__))
                    λ = (DynamicPPL.getconditioned_nested)(__context__, var"##vn#378")
                end
                (var"##value#380", __varinfo__) = (DynamicPPL.tilde_observe!!)(__context__, (DynamicPPL.check_tilde_rhs)(var"##dist#381"), λ, var"##vn#378", __varinfo__)
                var"##value#380"
            end
        end
        #= REPL[32]:3 =#
        begin
            var"##dist#386" = truncated(Normal(0, 1); lower = 0.0)
            var"##vn#383" = (DynamicPPL.resolve_varnames)((AbstractPPL.VarName){:error}(), var"##dist#386")
            var"##isassumption#384" = begin
                    if (DynamicPPL.contextual_isassumption)(__context__, var"##vn#383")
                        if !((DynamicPPL.inargnames)(var"##vn#383", __model__)) || (DynamicPPL.inmissings)(var"##vn#383", __model__)
                            true
                        else
                            error === missing
                        end
                    else
                        false
                    end
                end
            if (DynamicPPL.contextual_isfixed)(__context__, var"##vn#383")
                error = (DynamicPPL.getfixed_nested)(__context__, var"##vn#383")
            elseif var"##isassumption#384"
                begin
                    (var"##value#387", __varinfo__) = (DynamicPPL.tilde_assume!!)(__context__, (DynamicPPL.unwrap_right_vn)((DynamicPPL.check_tilde_rhs)(var"##dist#386"), var"##vn#383")..., __varinfo__)
                    error = var"##value#387"
                    var"##value#387"
                end
            else
                if !((DynamicPPL.inargnames)(var"##vn#383", __model__))
                    error = (DynamicPPL.getconditioned_nested)(__context__, var"##vn#383")
                end
                (var"##value#385", __varinfo__) = (DynamicPPL.tilde_observe!!)(__context__, (DynamicPPL.check_tilde_rhs)(var"##dist#386"), error, var"##vn#383", __varinfo__)
                var"##value#385"
            end
        end
        #= REPL[32]:4 =#
        y_hat = [mymodel(λ) for i = 1:100]
        #= REPL[32]:5 =#
        begin
            #= /home/raf/.julia/packages/DynamicPPL/txq74/src/compiler.jl:555 =#
            var"##retval#393" = begin
                    var"##dist#391" = MvNormal(y_hat, error)
                    var"##vn#388" = (DynamicPPL.resolve_varnames)((AbstractPPL.VarName){:y}(), var"##dist#391")
                    var"##isassumption#389" = begin
                            if (DynamicPPL.contextual_isassumption)(__context__, var"##vn#388")
                                if !((DynamicPPL.inargnames)(var"##vn#388", __model__)) || (DynamicPPL.inmissings)(var"##vn#388", __model__)
                                    true
                                else
                                    y === missing
                                end
                            else
                                false
                            end
                        end
                    if (DynamicPPL.contextual_isfixed)(__context__, var"##vn#388")
                        y = (DynamicPPL.getfixed_nested)(__context__, var"##vn#388")
                    elseif var"##isassumption#389"
                        begin
                            (var"##value#392", __varinfo__) = (DynamicPPL.tilde_assume!!)(__context__, (DynamicPPL.unwrap_right_vn)((DynamicPPL.check_tilde_rhs)(var"##dist#391"), var"##vn#388")..., __varinfo__)
                            y = var"##value#392"
                            var"##value#392"
                        end
                    else
                        if !((DynamicPPL.inargnames)(var"##vn#388", __model__))
                            y = (DynamicPPL.getconditioned_nested)(__context__, var"##vn#388")
                        end
                        (var"##value#390", __varinfo__) = (DynamicPPL.tilde_observe!!)(__context__, (DynamicPPL.check_tilde_rhs)(var"##dist#391"), y, var"##vn#388", __varinfo__)
                        var"##value#390"
                    end
                end
            #= /home/raf/.julia/packages/DynamicPPL/txq74/src/compiler.jl:556 =#
            @show y
            return (var"##retval#393", __varinfo__)
        end
    end
end
begin
    $(Expr(:meta, :doc))
    function test_model(y; )
        #= REPL[32]:1 =#
        return (DynamicPPL.Model)(test_model, NamedTuple{(:y,)}((y,)); )
    end
end
