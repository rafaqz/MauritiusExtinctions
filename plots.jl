
function plottimeline(speciesdata, island, category)
    timeline = getproperty(getproperty(pops, island), category)
    periods = timeline[!, 1]
    period_ranges = map(split.(periods, '-')) do (f, l)
        parse(Int, f):parse(Int, l)-1
    end
    period_intervals = map(split.(periods, '-')) do (f, l)
        Interval{:closed,:open}(parse(Int, f), parse(Int, l))
    end
    timerange = first(first(period_ranges)):last(last(period_ranges))
    groupnames = sort(unique(speciesdata.group))
    cschemes = (:purple, :dodgerblue, :orange, :green4, :yellow3, :firebrick)
    groupcolors = Dict(groupnames .=> cschemes[1:length(groupnames)])
    species_traits = filter(r -> !ismissing(r.LostLand_name) && r.LostLand_name in names(timeline), speciesdata)
    selected_timelines = timeline[!, in.(names(timeline), Ref(species_traits.LostLand_name))]
    order = [findfirst(r -> r.LostLand_name == n, Tables.rows(species_traits)) for n in names(selected_timelines) ]
    selected_species_traits = species_traits[order, :]
    abundance = zeros(Union{Int,Missing}, 
        Dim{:species}(names(selected_timelines)),
        Ti(timerange), 
    )
    for (i, interval) in enumerate(period_intervals)
        abundance[species=:, Ti=interval] .= Vector(selected_timelines[i, :])
    end
    groups = selected_species_traits.group
    groupheatmaps = map(groupnames) do gn
        kw = if gn == last(groupnames)
            (; xguide="") 
        else
            (; xguide="", xaxis=nothing) 
        end
        isingroup = groups .== gn
        scheme = cgrad([:lightgrey, groupcolors[gn]])
        any(isingroup) || return heatmap(; yticks=:none, showaxis=false, kw...)
        # TODO split up terrestrial and aquatic birds and animals
        ab = abundance[species=isingroup]
        perm = sortperm(selected_species_traits[isingroup, :], :Location)
        mask = falses(dims(abundance))
        if category == :native
            map(eachslice(ab, dims=:species), eachslice(mask, dims=:species)) do a, m
                any(x -> !ismissing(x), a) || return
                for i in eachindex(a)
                    !ismissing(a[i]) && break
                    a[i] = -1
                end
            end
        end
        p = heatmap(ab;
            c=scheme,
            clims=(-1, 3), colorbar=:none, 
            yguide=gn, yguidefontsize=14, yguide_position=:left,
            guidefontalign=:left, ytickdirection=:none, grid=:none,
            tickfontsize=9, ytickfontvalign=:vcenter,
            kw...
        )
    end
    animal_plot_fraction = 0.85
    animal_plot_heights = [count(==(g), groups) + 1 for g in groupnames] ./ 
        (length(groups) + length(groupnames)) .* animal_plot_fraction
    human_pop_plot = plot(getproperty(human_pop_timeline, island); 
        title=titlecase("$island $category"),
        c=:black, fill=(0, :grey),
        yguidefontrotation=-90,
        legend=:topleft, ticks=:none, xguide=""
    )
    sugar_plot = plot(sugar_timeline; 
        c=:navajowhite3, fill=(0, :beige),
        yguidefontrotation=-90, yguidefonthalign=:left,
        legend=:topleft, ticks=:none, xguide=""
    )
    other_plots = ()# (human_pop_plot, sugar_plot)
    other_plot_heights = map(_ -> (1 - animal_plot_fraction) / length(other_plots), other_plots)
    heights = [other_plot_heights..., animal_plot_heights...]
    plot(other_plots..., groupheatmaps...; 
        layout=grid(length(heights), 1; heights),
        size=(2000, 1000),
        link=:x,
    )
end

