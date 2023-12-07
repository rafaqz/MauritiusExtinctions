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

includet("optimisation.jl")
include("species_rules.jl")
include("raster_common.jl")

function gpu_cleanup(A)
    x, y, ti = lookup(A, (X, Y, Ti))
    Rasters.set(A, 
        Ti => Sampled(first(ti) - 0.5:last(ti) - 0.5; sampling=Intervals(Start())),
        X => LinRange(first(x), last(x), length(x)), 
        Y => LinRange(first(y), last(y), length(y)),
    )
end
uncleared = gpu_cleanup(Rasters.modify(BitArray, Raster("uncleared.nc")))
# forested = gpu_cleanup(modify(BitArray, Raster("forested.nc")))

pred_df = CSV.read("animals.csv", DataFrame)
introductions_df = CSV.read("introductions.csv", DataFrame)
mascarine_species_csv = "/home/raf/PhD/Mascarenes/Tables/mascarine_species.csv"
# @async run(`libreoffice $mascarine_species_csv`)
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
aggfactor = 8

lc_predictions_path = "$outputdir/lc_predictions.nc"
# netcdf has the annoying center locus for time
lc_predictions = RasterStack(lc_predictions_path) |>
    x -> rebuild(Rasters.modify(BitArray, x); missingval=false) |>
    x -> Rasters.set(x, Ti => Int.(maybeshiftlocus(Start(), dims(x, Ti), )))
aux = lc_predictions
include("species_rules.jl")
(; ruleset, pred_ruleset, inits, pred_inits, outputs, pred_outputs, outputs_kw, ag_masks, carrycap) = def_syms(
    pred_df, introductions_df, aggfactor, dems, masks, slope_stacks, island_extinct_tables, aux
);
# Rasters.rplot(aux[Ti=270])
island = :mus
@time sim!(outputs.mus, ruleset; proc=ThreadedCPU());
# @time sim!(pred_outputs.mus, pred_ruleset; proc=SingleCPU());
max_pops = map(outputs.mus) do frame
    map(keys(first(frame.pred_pops))) do key
        maximum(x -> getproperty(x, key), frame.pred_pops)
    end
end |> maximum

mkoutput = mk(inits[island], ruleset; carrycaps=max_pops, outputs_kw[island]...)
display(mkoutput)
Rasters.rplot(getproperty.(mkoutput[end].pred_pops, :mouse))
maximum(getproperty.(mkoutput[end].pred_pops, :mouse))
# vecmax(a, x) = max.(a, x)
# popmaxs = min.(60, mapreduce(s -> reduce(vecmax, s.pred_pops), vecmax, pred_outputs.mus))
masks.mus

using Stencils, StaticArrays, BenchmarkTools
rand(typeof(w)) .* rand(typeof(w))
Moore{2,2}(rand(24))
Moore(; radius=2)
A = StencilArray(rand(10, 10), Moore(2));
m = stencil(A, (2, 2))
mk_pred(pred_inits[island], pred_ruleset, ag_masks.mus; outputs_kw[island]...)
mk(inits[island], ruleset; outputs_kw[island]...)

# Random Spread
# approches zero as N approaches Inf
# but allows small populations to still spread because
# α increases with small N, so specific cells are more likely
# to be 1 or above
w = Window{2}(rand(5, 5))
N = 10
M = 3
α = (M / N)
A = rand(typeof(w)) .^ α 
scalars = (A .* (length(A) / sum(A)))
           .- 1) ./ N .+ 1
maximum(scalars)

length(pred_kernels[1].kernel)
x = zeros(7, 7)
kernel1 = DispersalKernel(
    stencil=moore,
    formulation=ExponentialKernel(100.0)
)
kernel2 = DispersalKernel(
    stencil=moore,
    formulation=ExponentialKernel(1.0)
)
kernel1
kernel1.kernel
kernel2.kernel
foreach(DynamicGrids.indices(kernel, (4, 4)), kernel.kernel) do i, k
    x[i...] = k
end
p = Makie.heatmap(x)
Colorbar(p.figure[1, 2], p.plot)


island_extinct_tables.mus.Hunting_preference

function getf()
    a = 1
    x -> begin
        a = 2
        a * x
    end
end
f = getf()
f.a

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
