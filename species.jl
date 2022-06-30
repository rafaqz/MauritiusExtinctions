using XLSX, DataFrames
using Plots, StatsPlots
using DimensionalData, DimensionalData.LookupArrays
using OrderedCollections
using IntervalSets

island_species = (
    mus=(native="Appendix 2", alien="Appendix 5"),
    reu=(native="Appendix 3", alien="Appendix 6"),
    rod=(native="Appendix 4", alien="Appendix 7"),
)

as_dataframe(xl, name::Symbol) = as_dataframe(xl, sheetnames[name])
function as_dataframe(xl, name::String)
    sheet = xl[name]
    df = DataFrame(XLSX.eachtablerow(sheet))
    df.Species .= strip.(isnumeric, df.Species)
    return df
end

# Turn classifications from the book into population values 0:3,
# representing extinct/not arrived, rare, common, abundant.
#
# Its not clear how to deal with grouped categories like `rats`

function filter_population(table)
    periods = names(table)[2:end-1]
    times = map(periods) do p
        s = split(p, '-')
        Base.parse(Int, s[1]), Base.parse(Int, s[2])
    end
    boundsmatrix = reinterpret(reshape, Int, times)
    timedim = Ti(Sampled(map(first, times);
        order=ForwardOrdered(),
        span=Explicit(boundsmatrix),
        sampling=Intervals(Start())),
    )
    populations = Dict{Symbol,Any}()
    for i in 1:size(table, 1)
        popvals = Vector{Union{Missing,Int}}(undef, length(periods))
        popvals .= missing
        started = false
        rawdata = table[i, 2:end-1]
        common_name = table[i, 1]
        local lastval = missing
        for j in eachindex(popvals)
            !started && ismissing(rawdata[j]) && continue
            rawval = rawdata[j]
            rawval = rawval isa String ? rawval[1:1] : rawval
            if haskey(data_key, rawval)
                category = data_key[rawval]
            else
                @warn "unidentified value in table: $rawval, for $common_name. `missing` used instead"
                continue
            end
            pop_estimate = if ismissing(category)
                lastval
            elseif category == "abundant" 
                3
            elseif category == "common" 
                2
            elseif category == "rare" 
                1
            elseif category == "uncertain" 
                lastval
            elseif category == "extinction/absence"
                0
            elseif category == "observed"
                @show common_name
                ismissing(lastval) ? 1 : lastval
            elseif category == "likely but unconfirmed"
                1
            elseif category == "several species unseparated"
                ismissing(lastval) ? 1 : lastval
            elseif category == "present no record"
                lastval
            elseif category == "not reported"
                lastval
            elseif category == "recorded introduction"
                1
            elseif category == "approximate introducion"
                1
            elseif category == "captive only"
                0
            elseif category == "I dont know what this is"
                # @info "category switch in $(table[i, 1])"
                # :movein
                missing
            elseif category == "move out of category" 
                # @info "category switch in $(table[i, 1])"
                # :moveout
                missing
            elseif category == "move into category"
                # @info "category switch in $(table[i, 1])"
                # :movein
                missing
            else
                continue
            end
            lastval = pop_estimate
            if ismissing(popvals[j])
                popvals[j] = pop_estimate
            end
            started = true
        end
        for i in reverse(eachindex(popvals))
            if ismissing(rawdata[i])
                popvals[i] = missing
            end
        end
        populations[Symbol(common_name)] = DimArray(popvals, timedim)
    end
    return DataFrame(:period => periods, populations...)
end

reverse_data_key = OrderedDict(
    "abundant" => "a",
    "common" => "b",
    "rare" => "c",
    "uncertain" => "d",
    "extinction/absence" => "e",
    "observed" => "f",
    "likely but unconfirmed" => "g",
    "several species unseparated" => "h",
    "present no record" => "i",
    "not reported" => "j",
    "recorded introduction" => "k",
    "approximate introduction" => "l",
    "captive only" => "m",
    "I dont know what this is" => "L",
    "move out of category"  =>  "N",
    "move into category" => "O",
    missing => missing,
)
data_key = OrderedDict(reverse(p) for p in reverse_data_key)

using CSV
species = CSV.File("/home/raf/PhD/Mauritius/mascarine_species.csv") |> DataFrame
lostland_names = collect(skipmissing(species[!, :LostLand_name]))
island_species = names(pops.mus.native)
map(x -> !in(x, lostland_names) ? x : missing, island_animals) |> skipmissing |> collect
# Load the transcribed XL file and turn each sheet into a dataframe
xlfile = "/home/raf/PhD/Mauritius/Data/LostLand/Mauritius_Lost Land of the Dodo_tables_translated symbols.xlsx"
# run(`libreoffice $xlfile`)
xl = XLSX.readxlsx(xlfile)
# filter(x -> x[1] == "Dodo", df)
pops = map(island_species) do island
    map(island) do sheetname
        @show sheetname
        filter_population(as_dataframe(xl, sheetname))
    end
end

timeline = pops.mus.alien
periods = timeline[!, 1]
groupnames = sort(unique(species.group))
cschemes = (:Purples, :Blues, :Oranges, :Greens, :Greys, :Reds)
groupcolors = Dict(groupnames .=> cschemes[1:length(groupnames)])
species_traits = filter(r -> !ismissing(r.LostLand_name) && r.LostLand_name in names(timeline), species)
select_timelines = timeline[!, in.(names(timeline), Ref(species_traits.LostLand_name))]
order = [findfirst(r -> r.LostLand_name == n, Tables.rows(species_traits)) for n in names(select_timelines) ]
selected_species_traits = species_traits[order, :]
abundance = Matrix(select_timelines)
# groupnums = [groupkey[g] for g in grouped.group]
# sorted_names = species_names[order]
# sorted_abundance = abundance[:, order]
groups = selected_species_traits.group
groupheatmaps = map(groupnames) do gn
    isingroup = groups .== gn
    labels = names(select_timelines)[isingroup]
    selected_abundance = abundance[:, isingroup]
    kw = gn == first(groupnames) ? (;) : (; yaxis=nothing) 
    heatmap(labels, periods, selected_abundance; 
        xguide=gn, c=groupcolors[gn], colorbar=:none,
        xrotation=45, xticks=(eachindex(labels), labels),
        tickfontsize=10, xtickfontvalign=:bottom, xtickdirection=:none,
        yflip=true, kw...
    )
end
widths = [count(==(g), groups) for g in groupnames] ./ length(groups)
plot(groupheatmaps...; 
    layout=grid(1, length(widths); widths),
    size=(2000, 2000),
)

species = map(pops) do island
    dfn = island.native
    extinction = map(names(dfn)) do name
        obs = dfn[!, name]
        i = findlast(x -> (!ismissing(x) && x != 0), obs)
        name => (isnothing(i) ? missing : dims(obs, Ti())[i])
    end |> OrderedDict
    dfa = island.alien
    introduction = map(names(dfa)) do name
        obs = dfa[!, name]
        i = findfirst(x -> !ismissing(x), obs)
        name => (isnothing(i) ? missing : dims(obs, Ti())[i])
    end |> OrderedDict
    (; extinction, introduction)
end


species_years = map(species) do island
    map(island) do class
        filter(sort(collect(class); by=last)) do x
            !ismissing(x[2]) && x[2] < 2000
        end
    end
end

declen = 10
decades = 1590:decade:1990
groups = map(species_years) do island
    map(island) do class
        groups = []
        for d in decades
            group = [s for s in class if s[2] in d..(d+declen-1)]
            push!(groups, group)
        end
        groups
    end
end

decade_counts = map(groups) do island
    map(class -> length.(class), island)
end
decade_annotations = map(groups) do island
    map(island) do class
        map(class) do x
            names = first.(x)
            join(names, ", ")
        end
    end
end
labels = [string.(keys(groups))...;;]
count_matrxs = map((introduction=:introduction, extinction=:extinction)) do key
    hcat(map(x -> getfield(x, key), values(decade_counts))...)
end

valigns = (mus=:bottom, reu=:center, rod=:top)
anns = map(decade_counts, decade_annotations, valigns) do cs, ts, valign
    offset = 2.5
    shift = valign == :top ? offset : valign == :bottom ? -offset : 0
    map(cs, ts) do c, t
        t = text.(t; pointsize=6, halign=:left)
        tuple.(0, decades .+ shift, t) 
    end
end
ann = anns.mus
kw = (linealpha=0, orientation=:horizontal, yflip=true, labels)
int = groupedbar(decades, count_matrxs.introduction; 
    xlab="Introductions", legend=:none, kw...
)
annotate!(int, anns.mus.introduction)
annotate!(int, anns.reu.introduction)
annotate!(int, anns.rod.introduction)
ext = groupedbar(decades, count_matrxs.extinction;
    xlab="Extinctions", yticks=:none, kw...
)
annotate!(ext, anns.mus.extinction)
annotate!(ext, anns.reu.extinction)
annotate!(ext, anns.rod.extinction)
plot(int, ext; size=(1000, 2000))#; layout=(2, 1))

# Rough ranking of commonness over time
ranked = map(pops) do origins
    map(origins) do origin
        sort(names(origin) .=> sum.(skipmissing.(eachcol(origin))); by=last, rev=true)
    end
end
ranked.reu.alien


key_invasives = [
    "goats",
    "cattle",
    "pigs",
    "cats",
    "rats",
    "Norway Rat",
    "Ship Rat",
    "Crab-eating Macacque",
]

# x = pops[:mauritius_invasive][!, "Black-naped Hare"]
pops.mus.alien[!, "cats"][At(1740)]
p = Plots.plot(Matrix(pops.rod.alien); labels=permutedims(names(pops.mus.alien)), size=(1200, 1200))

A = Matrix(pops.mus.alien[!, key_invasives])
p = Plots.plot(A .+ ((1:size(A, 2))./50)';
    labels=permutedims(key_invasives), opacity=0.5,
    legend=:bottomleft,
)

x = "Norway Rat"
Plots.plot(pops.mus.alien[!, x]; labels=x, size=(1000, 1000), ylims=(0, 3))



# Species distributions
#
# GBIF alternatives
# Arctos, neotoma, vertnet
using GBIF, Plots, Rasters, IntervalSets

# Google search all species, opening the next when the browser is closed.
for row in eachrow(subset(species, :Origin => ByRow(==("Endemic"))))
    ismissing(row["Common name"]) || ismissing(row["Species"]) && continue
    placenamed = false
    for place in ("RéuReunion", "Mauritius", "Rodrigues", "Mascarene")
        placenamed = occursin(place, row["Common name"]) 
        @show row["Common name"] place placenamed
        placenamed && break
    end
    placenamed && continue
    search = replace(row["Species"], " " => "+")
    run(`chromium https\://www.google.com/search\?q=$(search)`)
    # run(`chromium https\://www.google.com/search\?tbm=isch\&q=$(search)`)
end

for sp in species.Species
    ismissing(sp) && continue
    search = replace(sp, " " => "+")
    run(`chromium https\://www.google.com/search\?q=$(search)`)
    # run(`chromium https\://www.google.com/search\?tbm=isch\&q=$(search)`)
end

# Check in Hume
humepath = "/home/raf/zotero_library/Journal Article/undefined/2013/Hume_2013_A*"
for sp in species.Species
    ismissing(sp) && continue
    run(`evince --find=$(sp) "/home/raf/zotero_library/Journal Article/undefined/2013/Hume_2013_A synopsis of the pre-human avifauna of the Mascarene Islands.pdf"`)
end


endemics = DataFrames.subset(species, :Origin => ByRow(==("Endemic")); skipmissing=true)
lats, lons = (-22.0, -18.0), (55.0, 58.0)
bounds_flags = "decimalLatitude" => lats, "decimalLongitude" => lons

records = map(endemics[!, :Species]) do sp
    isnothing(sp) || ismissing(sp) && return sp => missing
    taxon = GBIF.taxon(sp)
    return sp => isnothing(taxon) ? missing : DataFrame(GBIF.occurrences(taxon, "limit"=>300, bounds_flags...))
end |> Dict

for (k, v) in records
    ismissing(v) && continue
    df = DataFrame(v)
    if all(ismissing.(df[!, :latitude]))
        records[k] = missing
    else
        @show k length(collect(skipmissing(df.latitude)))
    end
end

df = DataFrame(records["Gallinula chloropus"])
points = collect(zip(skipmissing(df.longitude), skipmissing(df.latitude)))

prec = Raster(WorldClim{Climate}, :prec; month=1, res="30s")
prec = Raster(CHELSA{Climate}, :prec; month=1, res)
mascarines_prec = prec[X = Interval(lons...), Y=Interval(lats...)]
Plots.plot(mascarines_prec)
scatter!(points; opacity=0.5, markershape=:circle)
obs = occurrences(, "limit"=>300)

while length(obs) < size(obs)
    @show length(obs) size(obs)
    occurrences!(obs)
end


# Tortoise
# Flat island/Garbriel island had 6000 tortoises. p 205
tortoise_carry_cap = 6000 / fi_area

# Macaque
# How large were historic populations?
# From Sussman and Tattersall
# 25000 - 35000
# Prefere secondary forest!
macaque_1986_pop = 30000
