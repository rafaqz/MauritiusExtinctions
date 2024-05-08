using TableView
using GLMakie, TerminalPager
GLMakie.activate!()

include("species_common.jl")


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
        # Remove missing points and years
        !ismissing(lon) && !ismissing(lat) && !ismissing(r.year) &&
        # Remove points not in mainland mauritis or reunion
        (extract(masks.mus, (lon, lat)).mask === true || extract(masks.reu, (lon, lat)).mask === true) &&
        # Has more than 2 decimal places
        !(round(lon * 100) / 100 ≈ lon || round(lat * 100) / 100 ≈ lat)
    end
    nrow(ocs) > 0 ? common => ocs : missing
end |> skipmissing |> Dict

function plot_ocs(common, ocs)
    ocs = filter(r -> !ismissing(r.year), ocs)
    points = collect(zip(ocs.decimalLongitude, ocs.decimalLatitude))
    points = map(p -> Float64.(p), points)
    fig = Figure()
    ax = Axis(fig[1, 1]; aspect=DataAspect(), title=common)
    Makie.heatmap!(ax, dems.mus)
    Makie.heatmap!(ax, native_veg.mus; colormap=:reds)
    Makie.heatmap!(ax, dems.reu)
    Makie.heatmap!(ax, native_veg.reu; colormap=:reds)
    # Makie.heatmap!(ax, dems.rod)
    s = Makie.scatter!(ax, points; 
        colormap=:thermal, 
        color=map(identity, ocs.year), 
        # colorrange=extrema(ocs.year)
    )
    Colorbar(fig[1, 2], s)
    return fig
end

common = "Mauritius Olive White-eye"
plot_ocs(common, occ_dfs[common])
occ_dfs
tbl = island_endemic_tables.mus
lc = lc_predictions.mus
mask = masks.mus
extract(lc[Ti=End] , [(57.367977, -20.434295)])

occurs_in = map(eachrow(tbl)) do species
    haskey(occ_dfs, species.Common_name) || return missing
    occ_df = occ_dfs[species.Common_name]
    occs = collect(zip(occ_df.decimalLongitude, occ_df.decimalLatitude))
    lcv = view(lc, Ti=End)
    timelines = map(oc -> Rasters.extract(lcv, oc; index=true), occs)
    in_native = map(timelines, occ_df.year) do (; index, native), y
        if ismissing(native) || ismissing(index) || ismissing(y) || (!mask[index[1:2]...])
            missing
        else
            lc[X(index[1]), Y(index[2]), Ti=Contains(min(2020, y))]
        end
    end
end;
occurs_in[1]

lc_fractions = map(occurs_in) do oin
    ismissing(oin) && return missing
    nonmissing = collect(skipmissing(oin))
    length(nonmissing) == 0 && return missing
    mean(x -> NamedVector(x), nonmissing)
end
nobs = map(occurs_in) do oin
    ismissing(oin) ? 0 : length(collect(skipmissing(oin)))
end

list = DataFrame(:Common_name=>tbl.Common_name, :lc=>lc_fractions, :nobs=>nobs)
df = filter(r -> !ismissing(r.lc), list)
expanded_list = map(df.Common_name, df.nobs, df.lc) do Common_name, nobs, lc
    (; Common_name, nobs, NamedTuple(lc)...)
end |> DataFrame


sp = "Grey Tomb Bat"
sort(occ_dfs[sp], :year) |> pager

for (k, df) in occ_dfs
    @show any(!=("PRESENT"), df.occurrenceStatus)
end

occ_dfs["Mauritius Olive White-eye"].year
occ_dfs["Mascarene Free-tailed Bat"]
sort(collect(keys(occ_dfs)))[1:20]

for (common, ocs) in occ_dfs 
    ocs = filter(r -> !ismissing(r.year), ocs)
    println(common)
    println()
    nrow(ocs) > 0 || continue
    fig = plot_ocs(common, ocs)
    makie = display(fig) 
    save("images/gbif/$common.png", fig)
    # while makie.window_open[] sleep(0.1) end
end

# for (common, ocs) in occ_dfs 
#     println(common)
#     ocs |> pager
# end
