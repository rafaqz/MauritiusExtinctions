includet("roads.jl")
includet("cost_distance.jl")
# includet("map_file_list.jl")
includet("ports.jl")
includet("slope.jl")
includet("raster_common.jl")

files = get_map_files()
slices = make_raster_slices(masks, lc_categories)
slices.mus.timelines.cleared

travel_origins = map(dems, ports) do d, p 
    o = map(_ -> Inf * u"hr", d)
    # foreach(dems, travel_origins) do d, oa
    # end
    for port in p.major
        o[map(Near, reverse(port))...] = 0.0u"hr"
    end
    for port in p.minor
        o[map(Near, reverse(port))...] = 1.0u"hr"
    end
    o
end
elevation = map(dem -> dem .* u"m", dems)
resistance = map(road_series) do island_road_series
    map(island_road_series) do roads
        broadcast(roads) do r
            if ismissing(r) 
                Inf
            elseif r
                relative_movement_speed[:dirt_road]
            else
                relative_movement_speed[:cleared_land]
            end
        end
    end
end
# plot(resistance.mus; clims=(0.0, 1.0))

# A history of woods and forests of Mauritius:
# 1 league = 4.83 kilometers per day in forest
# forst_walking_speed = 

# resistance = map(
#     (mus=slices.mus.timelines.cleared.cleared_1810, reu=slices.reu.timelines.cleared.cleared_1815),
#     (mus=way_masks.mus[:ways_1813], reu=way_masks.reu[:ways_1804])) do cleared, ways
#     broadcast(cleared, ways) do c, w
#         if w
#             relative_movement_speed[:dirt_road]
#         elseif c
#             relative_movement_speed[:cleared_land]
#         else
#             relative_movement_speed[:forest]
#         end
#     end
# end

@time travel_times = map(travel_origins, elevation, resistance) do o, e, r_series
    map(r_series) do r
        cost_distance(o, e, r; cellsize=step(lookup(dems.mus, X)) * 111.0u"km")
    end
end
Plots.plot(travel_times.mus)

# plot(travel_times.reu; clims=(0u"hr", 25u"hr"), legend=false, ticks=:none)
# plot(travel_times.mus; clims=(0u"hr", 12u"hr"), legend=false, ticks=:none)

# anims = map(travel_times) do tt
#     @animate for (A, t) in zip(tt, lookup(tt, Ti))
#         plot(A; legend=false, title=string("travel time ", t), clims=(0, 25))
#     end
# end
# gif(anims.mus, "mus_travel_times.gif", fps=0.8)
# gif(anims.reu, "reu_travel_times.gif", fps=0.8)

# @time travel_times = map(travel_origins, elevation, resistance) do o, e, r
#     cost_distance(o, e, r; cellsize=step(lookup(dems.mus, X)) * 111.0u"km")
# end
# plot(map(t -> plot(t; clims=(0, 40)), travel_times)...)

# plot(plot(travel_times.reu), plot(slices.reu.timelines.cleared.cleared_1815), Plots.contourf(dems.reu; levels=[0:500:3000...]))
# reu_dist_cleared = replace_missing(travel_times.reu) .* rebuild(replace(slices.reu.timelines.cleared.cleared_1815, 0 => missing); missingval=missing)
# mus_dist_cleared = replace_missing(travel_times.mus) .* rebuild(replace(slices.mus.timelines.cleared.cleared_1810, 0 => missing); missingval=missing)
# mean(skipmissing(mus_dist_cleared))
# mean(skipmissing(reu_dist_cleared))
# plot(mus_dist_cleared; clims=(0, 10))
# plot(reu_dist_cleared; clims=(0, 10))
# count(_ -> true, skipmissing(mus_dist_cleared))
# count(_ -> true, skipmissing(reu_dist_cleared))
# size(mus_dist_cleared)
# size(reu_dist_cleared)
# plot(slope.mus)
