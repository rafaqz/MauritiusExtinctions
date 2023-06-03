using CSV, DataFrames, DimensionalData

includet("common.jl")

function interpolate_years(data, keys)
    pop_vec = map(all_years) do year
        years, pops = map(k -> getproperty(data, k), keys)
        i = searchsortedlast(years, year)
        if i == 0
            return 0.0
        elseif i >= length(pops)
            return Float64(last(pops))
        elseif years[i] == year
            return Float64(pops[i])
        else
            frac = (year - years[i]) / (years[i + 1] - years[i])
            return pops[i] + (pops[i + 1] - pops[i]) * frac
        end
    end
end

# Population
# human_pop = CSV.File(joinpath(workdir, "Data/Population/Population.csv")) |> DataFrame
sugar_cane = CSV.File(joinpath(datadir, "Population/Sugarcane.csv")) |> DataFrame

all_years = 1600:2020
mus_census_pop = CSV.File(joinpath(datadir, "Population/lutz_census.csv")) |> DataFrame
mus_early_pop = CSV.File(joinpath(datadir, "Population/early_population.csv")) |> DataFrame
mus_pop = reduce(vcat, [mus_early_pop, mus_census_pop])
reu_pop = CSV.File(joinpath(datadir, "Population/reunion_population.csv")) |> DataFrame

human_pop_timeline = map((mus=mus_pop, reu=reu_pop, rod=reu_pop)) do pop
    pop_vec = interpolate_years(pop, (:Year, :Population))
    DimArray(pop_vec, Ti(all_years); name=:Human_Population)
end
human_pop_timeline.reu[At(1700)]

sugar_vec = interpolate_years(sugar_cane, (:Year, :Area))
sugar_timeline = DimArray(sugar_vec, Ti(all_years); name=:Area_Sugar)

# Plot
# Plots.plot(human_pop_timeline.mus)
# Plots.plot(sugar_timeline)
