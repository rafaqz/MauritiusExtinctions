using XLSX, DataFrames
using DimensionalData, DimensionalData.LookupArrays

const sheetnames = (
    mauritius_native="Apendix 2", # Yes Apendix is misspelled
    reunion_native="Appendix 3",
    rogriguez_native="Appendix 4 ", # Space at end needed
    mauritius_invasive="Appendix 5",
    reunion_invasive="Appendix 6",
    rogrigues_invasive="Appendix 7",
)

as_dataframe(xl, name::String) = as_dataframe(xl, sheetnames[name])
function as_dataframe(xl, name::String)
    sheet = xl[name]
    return DataFrame(XLSX.eachtablerow(sheet))
end

reverse_data_key = Dict(
    "abundant" => "a",
    "common" => "b",
    "rare" => "c",
    "uncertain" => "d",
    "extinction" => "e",
    "observed" => "f",
    "unconfirmed" => "g",
    "several species unseparated" => "h",
    "present no record" => "i",
    "not reported" => "j",
    "recorded" => "k",
    "introduced" => "l",
    "captive only" => "m",
    "I dont know what this is" => "L",
    "move out of category"  =>  "N",
    "move into category" => "O",
    missing => missing,
)

data_key = Dict(reverse(p) for p in reverse_data_key)

function filter_population(table)
    periods = names(table)[3:end-1]
    times = map(periods) do p
        s = split(p, '-')
        parse(Int, s[1]), parse(Int, s[2])
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
        local lastval = missing
        for j in eachindex(popvals)
            rawdata = table[i, 3:end-1] 
            !started && ismissing(rawdata[j]) && continue
            # @show i j
            rawval = rawdata[j]
            rawval = rawval isa String ? rawval[1:1] : rawval
            if haskey(data_key, rawval)
                category = data_key[rawval]
            else
                @warn "unidentified value in table: $rawval. `missing` used instead"
                continue
            end
            pop_estimate = if category == "abundant" 
                1
            elseif category == "common" 
                2
            elseif category == "rare" 
                3
            elseif category == "extinction"
                0
            elseif category == "observed"
                ismissing(lastval) ? 1 : lastval
            elseif category == "unconfirmed"
                ismissing(lastval) ? 1 : lastval
            elseif category == "several species unseparated"
                ismissing(lastval) ? 1 : lastval
            elseif category == "present no record"
                ismissing(lastval) ? 1 : lastval
            elseif category == "not reported"
                lastval
            elseif category == "recorded"
                1
            elseif category == "introduced"
                1
            elseif category == "captive only"
                0
            elseif category == "I dont know what this is"
                @info "category switch in $(table[i, 1])"
                1
            elseif category == "move out of category" 
                @info "category switch in $(table[i, 1])"
                1
            elseif category == "move into category"
                @info "category switch in $(table[i, 1])"
            elseif category === missing
                missing
            end
            popvals[j] == pop_estimate
        end
        populations[table[i, 1]] = DimArray(popvals, timedim)
    end
    return populations
end



xlfile = "/home/raf/PhD/Mauritius/LostLand/Mauritius_Lost Land of the Dodo_tables_translated symbols.xlsx"
xl = XLSX.readxlsx(xlfile)
pops = map(sheetnames) do sheetname
    @show sheetname
    filter_population(as_dataframe(xl, sheetname))
end;
nothing

# using TableView, Blink
# using ProfileView
# @profview table = as_dataframe(xl, :mauritius_invasive)
# @profview w = Blink.Window()
# body!(w, TableView.showtable(table))

# pops[:mauritius_native]
