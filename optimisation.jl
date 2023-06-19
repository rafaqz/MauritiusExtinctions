using Statistics
# Simple non-spatial loss functions

# function loss_population(output, species_timeseries; carry_caps=nothing, npopclasses=3)
#     population_sums = map(o -> sum(o[:native], output))
#     population_classes = map(population_sums) do popsums
#         # Convert each population into abundance clases
#         # e.g. 0 = extinct, 1 = rare, 2 = common, 3 = abundant
#         ceil.(Int, popsums ./ carry_caps .* npopclasses), 
#     end
#     # Comparing population size classes directly is a simple way to compute loss.
#     # But probably we dont care if the population is common or abundant?
#     # Rare and extinct are the interesting parts
#     # Also it would be better for the optimiser to maintain a
#     # smooth gradient, rather than discretizing population like this.
#     loss = sum(count.(==, population_adjusted, species_timeseries))
    
#     # TODO: spatial components?
#     return loss
# end

# Multiple sims - get mean counts
function presence_loss(predictions::Vector{<:Vector{<:StaticVector}}; kw...)
    map(predictions) do p
        presence_loss(p; kw...)
    end |> mean
end
# Single sim
function presence_loss(predictions::Vector{<:StaticVector}; last_obs, tspan, ncells)
    losses = map(tspan, predictions) do t, ps
        map(ps, last_obs) do p, l_obs
            if t <= l_obs 
                # Not extinct yet so there should still be some presences, 
                # penalise time before last observation one per year
                p == zero(p) ? (l_obs - t + 1.0) : 0.0
            else
                # Penalise presences and time since last observation
                p > zero(p) ? p * (t - l_obs) : 0.0
            end |> Float64
        end
    end

    return sum(losses)
end
