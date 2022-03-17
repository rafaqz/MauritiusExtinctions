using XLSX, DataFrames
using Plots
using DimensionalData, DimensionalData.LookupArrays

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
    periods = names(table)[3:end-1]
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
    populations = Dict()
    for i in 1:size(table, 1)
        popvals = Vector{Union{Missing,Int}}(undef, length(periods))
        popvals .= missing
        started = false
        rawdata = table[i, 3:end-1] 
        local lastval = missing
        for j in eachindex(popvals)
            !started && ismissing(rawdata[j]) && continue
            rawval = rawdata[j]
            rawval = rawval isa String ? rawval[1:1] : rawval
            if haskey(data_key, rawval)
                category = data_key[rawval]
            else
                @warn "unidentified value in table: $rawval, for $(table[i, 1]). `missing` used instead"
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
                ismissing(lastval) ? missing : lastval
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
        name = strip(isnumeric, table[i, 1])
        populations[name] = DimArray(popvals, timedim)
    end
    return DataFrame(populations)
end

reverse_data_key = Dict(
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

data_key = Dict(reverse(p) for p in reverse_data_key)

# Load the transcribed XL file and turn each sheet into a dataframe
xlfile = "/home/raf/PhD/Mauritius/Data/LostLand/Mauritius_Lost Land of the Dodo_tables_translated symbols.xlsx"
# run(`libreoffice $xlfile`)
xl = XLSX.readxlsx(xlfile)
# df = as_dataframe(xl, :ms_native)
# df = as_dataframe(xl, :mus_invasive)
# filter(x -> x[1] == "Dodo", df)

pops = map(island_species) do island
    map(island) do sheetname
        @show sheetname
        filter_population(as_dataframe(xl, sheetname))
    end
end
ranked = map(pops) do origins
    map(origins) do origin
        sort(names(origin) .=> sum.(skipmissing.(eachcol(origin))); by=last, rev=true)
    end
end
ranked.reu.alien

key_invasives = [
    # "goats",
    # "pigs",
    # "cats",
    "rats",
    "Norway Rat",
    "Ship Rat",
    # "Crab-eating Macacque",
]

# x = pops[:mauritius_invasive][!, "Black-naped Hare"]
pops.mus.alien[!, "cats"][At(1740)]

p = Plots.plot(Matrix(pops.rod.alien); labels=permutedims(names(pops.mus.alien)), size=(1200, 1200))
A = Matrix(pops.mus.alien[!, key_invasives])
p = Plots.plot(A .+ rand(size(A)...) .* 0.1;
    labels=permutedims(key_invasives), opacity=0.5,
)
x = "Norway Rat"
Plots.plot(pops.mus.alien[!, x]; labels=x, size=(1000, 1000))



# Species distributions
#
# GBIF alternatives
# Arctos, neotoma, vertnet
using GBIF, Plots, Rasters, IntervalSets
species = CSV.File("/home/raf/PhD/Mauritius/mascarine_species.csv") |> DataFrame

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
