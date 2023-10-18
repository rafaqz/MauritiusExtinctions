using ModelParameters
using StaticArrays
using DynamicGrids
using Dispersal
using LandscapeChange
using CSV
using DataFrames
using XLSX
using TerminalPager
using Rasters
using GLMakie
using NCDatasets
using Geomorphometry

includet("/home/raf/.julia/dev/LandscapeChange/src/makie.jl")
includet("optimisation.jl")
include("species_rules.jl")
include("raster_common.jl")

uncleared_raw = modify(BitArray, Raster("uncleared.nc"))
x, y = lookup(uncleared_raw, (X, Y))
uncleared = set(uncleared_raw, 
    Ti => 1600:2018, 
    X => LinRange(first(x), last(x), length(x)), 
    Y => LinRange(first(y), last(y), length(y))
)

pred_df = CSV.read("animals.csv", DataFrame)
mascarine_species_csv = "/home/raf/PhD/Mascarenes/Tables/mascarine_species.csv"
@async run(`libreoffice $mascarine_species_csv`)
all_species = CSV.read(mascarine_species_csv, DataFrame) |> 
    x -> subset(x, :Species => ByRow(!ismissing))
island_tables = map((mus=:mus, reu=:reu, rod=:rod)) do key
    df = DataFrame(subset(all_species, key => x -> .!ismissing.(x)))
    df.presence = df[!, key]
    df.extinct = map(df[!, "$(key)_extinct"]) do e
        ismissing(e) ? missing : eval(Meta.parse(e))::UnitRange
    end
    df.introduced = map(df[!, "$(key)_introduced"]) do e
        ismissing(e) ? missing : eval(Meta.parse(e))::UnitRange
    end
    df
end
island_extinct_tables = map(island_tables) do tbl
    DataFrame(subset(tbl, :extinct => x -> .!ismissing.(x)))
end

get_species_names(table) = Tuple(Symbol.(replace.(skipmissing(table.Species), Ref(" " => "_"))))
island_names = NamedTuple{keys(island_tables)}(keys(island_tables))
island_extinct = map(island_tables, island_names) do table, name
    subset(table, Symbol(name, "_extinct") => ByRow(!ismissing))
end
island_extinct_names = map(get_species_names, island_extinct)
include("species_rules.jl")
aggfactor = 8
ruleset, inits, outputs, outputs_kw = def_syms(
    pred_df, aggfactor, dems, masks, slope_stacks, island_extinct_tables
);

island_extinct_tables.Hunting_preference

using ProfileView
using BenchmarkTools
using CUDA
using Adapt
CUDA.allowscalar(false)
CUDA.allowscalar(true)
@profview 
@time sim!(outputs.mus, ruleset; proc=SingleCPU());
sum(outputs.mus[end].pred_pops)
@time sim!(outputs.mus, ruleset; proc=CPUGPU());
# @time sim!(outputs.mus, ruleset; proc=ThreadedCPU(), opt=SparseOpt());
@time for i in 1:10
    sim!(outputs.mus, ruleset; proc=CPUGPU());
end
cu_output = Adapt.adapt(CuArray, outputs.mus)
@time sim!(cu_output, ruleset; proc=CuGPU());

function f(a, b) 
    Threads.@threads for i in 1:10000000000
        a[] = b[] + 1
        b[] = b[] + a[]
    end
    return b
end
@time f(Ref(1), Ref(2))

pairs(outputs.mus[213])
outputs.mus[At(1800)] |> pairs

island = :mus
mk(inits[island], ruleset; outputs_kw[island]...)


@time map(outputs) do output
    sim!(output, ruleset);
end

# @time sims = let trans_output=trans_output, ruleset=ruleset 
#     tasks = map(1:100) do i
#         Threads.@spawn sim!(deepcopy(trans_output), ruleset).frames
#     end
#     map(fetch, tasks)
# end;

# l = presence_loss(sims; tspan, ncells, last_obs)
# pairs(l)


# using CUDA, Adapt
# CUDA.allowscalar(false)
# using ProfileView @profview 
# cu_output = Adapt.adapt(CuArray, output);
# try
    # sim!(cu_output, ruleset; proc=CuGPU());
# catch err
    # code_typed(err; interactive=true)
# end
#
