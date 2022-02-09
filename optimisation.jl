# Simple non-spatial loss functions


function loss_population(output, species_timeseries; carry_caps=nothing, npopclasses=3)
    population_sums = map(o -> sum(o[:native], output))
    population_classes = map(population_sums) do popsums
        # Convert each population into abundance clases
        # e.g. 0 = extinct, 1 = rare, 2 = common, 3 = abundant
        ceil.(Int, popsums ./ carry_caps .* npopclasses), 
    end
    # Comparing population size classes directly is a simple way to compute loss.
    # But probably we dont care if the population is common or abundant?
    # Rare and extinct are the interesting parts
    # Also it would be better for the optimiser to maintain a
    # smooth gradient, rather than discretizing population like this.
    loss = sum(count.(==, population_adjusted, species_timeseries))
    
    # TODO: spatial components?
    return loss
end

function loss_presence(output::Vector{Matrix{SVector{Bool}}}, observations::Vector{Vector{Bool}}; grid=:native_species)
    # The number of predictions is the number of species by the number of timesteps
    npredictions = length(observations) * length(observations[1])
    # Convert simulation grid to presence/absense predictions
    predictions = map(o[grid]) do m
        # Reinterpret each Matrix of SVector tuples as a 3d Arrays.
        A = reinterpret(reshape, eltype(A[1][1]))
        # Check if there are any prescences in the array slice for each species
        any(A; dims=(2, 3))::Vector{Bool}
    end::Vector{Vector{Bool}}
    # Reduce all species at all timesteps to a sum correct predictions
    ncorrect = mapreduce(+, predictions, observations) do species_predictions, species_observations
        # observations and predictions are Booleans, so we just & them together
        mapreduce(&, +, species_predictions, species_observations)
    end
    # The loss is the fraction of correct predictions
    loss = ncorrect / npredictions

    # TODO: spatial components?
    return loss
end

function loss_probability(output, species_timeseries, carry_caps)
    # Work in progress...
    population_sums = map(o -> sum(o[:native], output))
    population_classes = map(population_sums) do popsums
        # Convert each population into abundance clases
        # e.g. 0 = extinct, 1 = rare, 2 = common, 3 = abundant
        popsums ./ carry_caps .* npopclasses 
    end

    # TODO: spatial components?
    # For probabilities we can keep things in floating point?
    sum(count.(==, population_adjusted, species_timeseries))
end


A = [(true, false, true) for i in 1:10, j in 1:5]

re = reinterpret(reshape, eltype(A[1][1]), A)
any(re; dims=(2, 3))
