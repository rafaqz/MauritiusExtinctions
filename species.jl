using XLSX, DataFrames, TableView, Blink
using DimensionalData, DimensionalData.LookupArrays

xlfile = "/home/raf/PhD/Mauritius/LostLand/Mauritius_Lost Land of the Dodo_tables_translated symbols.xlsx"

xl = XLSX.readxlsx(xl)

sheetnames = (
    mauritius_native="Apendix 2", # Yes this is misspelled
    reunion_native="Appendix 3",
    rogriguez_native="Appendix 4 ",
    mauritius_invasive="Appendix 5",
    reunion_invasive="Appendix 6",
    rogrigues_invasive="Appendix 7",
)

as_dataframe(xl, name::String) = as_dataframe(xl, sheetnames[name])
function as_dataframe(xl, name::String)
    sheet = xl[name]
    return DataFrame(XLSX.eachtablerow(sheet))
end


w = Blink.Window()
body!(w, TableView.showtable(table))

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

pop_key = Dict(
    "a" => 1,
    "b" => 2,
    "c" => 3,
)

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
        popvals = [map(x -> get(pop_key, x, missing), table[2, 3:end-1])...]
        populations[table[i, 1]] = DimArray(popvals, timedim)
    end
    return populations
end

pops = map(sheetnames) do sheetname
    filter_population(as_dataframe(xl, sheetname))
end

pops[:mauritius_native]
