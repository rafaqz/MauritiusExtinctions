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

# includet("optimisation.jl")
include("species_rules.jl")
include("raster_common.jl")
include("makie.jl")
fig = Figure()

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

pred_df = CSV.read("tables/animals.csv", DataFrame)
introductions_df = CSV.read("tables/introductions.csv", DataFrame)
mascarene_species_csv = "tables/mascarene_species.csv"
# @async run(`libreoffice $mascarene_species_csv`)
all_species = CSV.read(mascarene_species_csv, DataFrame) |> 
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

rast = Raster(rand(10, 10, 5), (X(), Y(), Band([:a, :b, :c, :d, :e])))
stack = RasterStack(rast; layersfrom=Band)

lc_predictions_paths = (
    mus="$outputdir/lc_predictions_mus.nc",
    reu="$outputdir/lc_predictions_reu.nc",
    rod="$outputdir/lc_predictions_rod.nc",
)

# netcdf has the annoying center locus for time
lc_predictions = map(lc_predictions_paths) do path
    lc_predictions = RasterStack(path) |>
        x -> rebuild(Rasters.modify(BitArray, x); missingval=false) |>
        x -> Rasters.set(x, Ti => Int.(maybeshiftlocus(Start(), dims(x, Ti), )))
end
Rasters.rplot(lc_predictions.mus[Ti=At(2010)])

k = :reu
k = :rod
k = :mus
include("species_rules.jl")
(; ruleset, islands) = def_syms(pred_df, introductions_df, aggfactor, dems, masks, slope_stacks, island_extinct_tables, lc_predictions);
(; output, init, output_kw) = islands[k]
@time sim!(output, ruleset; proc=ThreadedCPU(), printframe=true);
Makie.plot(mean(map(i -> getproperty.(output[end-i].pred_pop, :macaque), 0:20)); colormap=:magma)
maximum(getproperty.(mkoutput[end].pred_pops, :mouse))

# Get the max color for Makie
max_pops = map(output) do frame
    map(propertynames(first(frame.pred_pop))) do key
        maximum(x -> getproperty(x, key), frame.pred_pop)
    end
end |> maximum

mkoutput = mk(init, ruleset; carrycaps=max_pops, output_kw...)
display(mkoutput)


function predict_extinctions(ruleset, init; tspan, kw...)
    output = TransformedOutput(init; tspan, kw...) do data
        map(>(0), sum(data.species))
    end
    sim!(output, ruleset)
    # Return NamedVector of years to extinction
    return sum(sim) .+ first(tspan) .- 1
end

function predict_extinctions(ruleset, inits; kw...)
    map(islands) do 
        predict_extinctions(ruleset, island.init; 
            tspan=island.tspan, aux=island.aux, kw...
        )
    end
end
