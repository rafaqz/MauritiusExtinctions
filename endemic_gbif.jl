include("species_tables.jl")
include("raster_common.jl")
using GLMakie, TerminalPager

# Retreive all data from GBIF
# using GBIF2
# for name in endemic_species.GBIFSpecies
#     println(name)
#     csvpath = "/home/raf/PhD/Mascarenes/MauritiusExtinctions/tables/occurrence/$name.csv"
#     isfile(csvpath) && continue
#     species = GBIF2.species_match(name)
#     isnothing(species) && continue
#     occurrences = occurrence_search(species; limit=10000)
#     CSV.write(csvpath, occurrences)
# end

occ_dfs = map(endemic_species.GBIFSpecies, endemic_species.Common_name) do name, common
    csvpath = "/home/raf/PhD/Mascarenes/MauritiusExtinctions/tables/occurrence/$name.csv"
    isfile(csvpath) || return missing
    df = CSV.read(csvpath, DataFrame)
    nrow(df) > 0 || return missing
    ocs = filter(df) do r
        lat = r.decimalLatitude
        lon = r.decimalLongitude
        !ismissing(lon) && !ismissing(lat) && lon in 50 .. 60 && lat in -25 .. -15
    end
    nrow(ocs) > 0 ? common => ocs : missing
end |> skipmissing |> Dict

function plot_ocs(common, ocs)
    ocs = filter(r -> !ismissing(r.year), ocs)
    points = collect(zip(ocs.decimalLongitude, ocs.decimalLatitude))
    points = map(p -> Float64.(p), points)
    fig = Figure()
    ax = Axis(fig[1, 1]; aspect=DataAspect(), title=common)
    p = Makie.heatmap!(ax, dems.mus)
    Makie.heatmap!(ax, dems.reu)
    # Makie.heatmap!(ax, dems.rod)
    s = Makie.scatter!(ax, points; 
        colormap=:thermal, 
        color=map(identity, ocs.year), 
        # colorrange=extrema(ocs.year)
    )
    Colorbar(fig[1, 2], s)
    return fig
end

occ_dfs["Mauritius Olive White-eye"].year
occ_dfs["Mauritius Olive White-eye"] |> plot_ocs
sort(collect(keys(occ_dfs)))


for (common, ocs) in occ_dfs 
    ocs = filter(r -> !ismissing(r.year), ocs)
    println(common)
    println()
    nrow(ocs) > 0 || continue
    fig = plot_ocs(common, ocs)
    makie = display(fig) 
    while makie.window_open[] 
        sleep(0.1)
    end
end

for (common, ocs) in occ_dfs 
    println(common)
    ocs |> pager
end
