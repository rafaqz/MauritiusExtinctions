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

as_dataframe(xl, name::Symbol) = as_dataframe(xl, sheetnames[name])
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
        popvals = Vector{Union{Missing,Symbol,Int}}(undef, length(periods))
        popvals .= missing
        started = false
        rawdata = table[i, 3:end-1] 
        local lastval = missing
        extinct = false
        for j in eachindex(popvals)
            !started && ismissing(rawdata[j]) && continue
            # @show i j
            rawval = rawdata[j]
            rawval = rawval isa String ? rawval[1:1] : rawval
            if haskey(data_key, rawval)
                category = data_key[rawval]
            else
                @warn "unidentified value in table: $rawval, for $(table[i, 1]). `missing` used instead"
                continue
            end
            pop_estimate = if extinct 
                0
            elseif ismissing(category)
                lastval
            elseif category == "abundant" 
                3
            elseif category == "common" 
                2
            elseif category == "rare" 
                1
            elseif category == "uncertain" 
                lastval
            elseif category == "extinction"
                extinct = true
                0
            elseif category == "observed"
                ismissing(lastval) ? :observed : lastval
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
                :movein
            elseif category == "move out of category" 
                # @info "category switch in $(table[i, 1])"
                :moveout
            elseif category == "move into category"
                # @info "category switch in $(table[i, 1])"
                :movein
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

# Load the transcribed XL file and turn each sheet into a dataframe
xlfile = "/home/raf/PhD/Mauritius/Data/LostLand/Mauritius_Lost Land of the Dodo_tables_translated symbols.xlsx"
xl = XLSX.readxlsx(xlfile)
pops = map(sheetnames) do sheetname
    @show sheetname
    filter_population(as_dataframe(xl, sheetname))
end;
df = as_dataframe(xl, :mauritius_native)
filter(x -> x[1] == "Dodo", df)


@info "TODO: how to classify extinction point when there is uncertainty"
x = pops[:mauritius_native][!, "Dodo"]
# using Plots
# plot(x)

# pops[:mauritius_invasive]

# using TableView, Blink
# using ProfileView
# @profview table = as_dataframe(xl, :mauritius_invasive)
# @profview w = Blink.Window()
# body!(w, TableView.showtable(table))

