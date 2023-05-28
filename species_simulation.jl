using DynamicGrids, Dispersal, Shapefile
using LandscapeChange
using CSV, DataFrames
using ModelParameters

includet("raster_common.jl")

# settlement = (
#     1598, "Port de Warwick/Grand port used as stopover", "Dutch", -20.3754, 57.7228,
#     1606, "Port Louie used as stopover", "Dutch", -20.1597, 57.5012,
#     1635, "First garrison at Port de Warwick", "Dutch", -20.3754, 57.7228, "Ebony clearing in this time",
#     1645, "Road build in Flacq",
#     1655, "Three settlements: Grand port bay, Flacq, Trou dEau Douce",
#     ~1700, "Flacq and Black River settlements",
#     1721, "Colony at port louie founded", "French", -20.1597, 57.5012,
#     1806, "Mahebourge founded", "French", -20.4045, 57.7028,
#     1814, "Mauritius ceded to britain", "English",
# )

s = CSV.read("/home/raf/PhD/Mascarenes/Tables/mascarine_species.csv", DataFrame) |> 
    x -> subset(x, :Species => ByRow(!ismissing))

island_tables = map((mus=:mus, reu=:reu, rod=:rod)) do key
    subset(s, key => x -> .!ismissing.(x))
end
map(length ∘ collect ∘ skipmissing, (s.mus_extinct, s.reu_extinct, s.rod_extinct))
# filter(e) do r
    # ismissing(r.rod) && !ismissing(r.rod_extinct)
# end

get_species_names(table) = Tuple(Symbol.(replace.(skipmissing(table.Species), Ref(" " => "_"))))
island_extinct = map(island_tables, island_names) do table, name
    subset(table, Symbol(name, "_extinct") => ByRow(!ismissing))
end
island_names = NamedTuple{keys(island_species)}(keys(island_species))
island_extinct_names = map(get_species_names, island_extinct)
species_params = map(s.Rmax, s.Max_Density, s.mus_extinct, s.reu_extinct, s.rod_extinct) do rmax, max_density, mus_extinct, reu_extinct, rod_extinct
    (; rmax, max_density, mus_extinct, reu_extinct, rod_extinct, hunting_suscept=1.0, cat_suscept=1.0, rat_suscept=1.0)
end |> Tuple |> NamedVector{get_species_names(s)}

island_params = map(island_extinct_names, island_names) do keys, island
    spec = species_params[keys]
    map(spec) do s
        extinct = s[Symbol(island, "_extinct")]
        base = (; s[(:rmax, :max_density, :hunting_suscept, :cat_suscept, :rat_suscept)]..., extinct)
    end
end

island_populations = map(island_params) do params
    v = NamedVector(map(_ -> 100.0f0, params))
    fill(v, 100, 100)
end

island_columns = map(island_params) do params
    k = keys(first(params))
    map(k) do key
        NamedVector(map(r -> r[key], params))
    end |> NamedTuple{k}
end

island_columns.mus.max_density
island_params.mus

init_pops = map(island_columns, masks) do params, mask
    fill(map(Float32, params.max_density), size(mask)) .* mask
end

tspan = 1600:2020

growth_rules = map(island_columns) do island
    growth = Dispersal.LogisticGrowth{:populations,:populations}(;
        rate=island.rmax,
        carrycap=island.max_density,
        timestep=1,
        nsteps_type=Float32,
    )
end

hunting = let hunting_suscept = hunting_suscept
    Cell{:populations,:populations}() do data, pops, I
        hp = get(data, Aux{:hunting_pressure}(), I)
        pops .* (1 .- hunting_suscept .* hunting_pressure)
    end
end

cat_predation = let cat_suscept = cat_suscept, 
                    weights = weights
    Cell{Tuple{:populations,:cats}}() do data, (pops, cats), I
        max_predation = pops .* cat_suscept .* cats
        max_nutrition = weights .* meat_kj .* potential_predation
        max_required = cats * cat_K .* max_nutrition
        predation = max_predation .* frac
        newpops .- predation, newpreds
    end
end

ruleset = Ruleset(growth, spread, cat_predation, rat_predation, hunting)

output = MakieOutput(init_state;
    aux,
    fps=10,
    tspan,
    store=false,
    mask=revmasks.mus,
    boundary=Remove(),
    padval=0,
    ruleset,
    sim_kw=(; printframe=true),
) do fig, frame
    axis = Axis(fig[1, 1])
    species = Observable(Array(frame[].landcover))
    i = Observable(1)
    on(frame) do f
        species[] = getindex.(f, i[])
        notify(landcover)
    end
    hm = Makie.image!(axis, landcover; colorrange=(-0.5, 5.5), colormap=:inferno)
    return nothing
end
