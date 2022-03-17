using CSV
using DataFrames
using Plots

# Population
# human_pop = CSV.File(joinpath(workdir, "Data/Population/Population.csv")) |> DataFrame
# sugar_cane = CSV.File(joinpath(workdir, "Data/Population/Sugarcane.csv")) |> DataFrame

census_pop = CSV.File(joinpath(datadir, "Population/lutz_census.csv")) |> DataFrame
early_pop = CSV.File(joinpath(datadir, "Population/early_population.csv")) |> DataFrame
combined_pop = reduce(vcat, [early_human_pop, census_human_pop])
all_years = 1600:2020
pop_vec = map(all_years) do year
    years, pops = human_population.Year, human_population.Population
    i = searchsortedlast(years, year)
    if years[i] == year
        return Float64(pops[i])
    else
        frac = (year - years[i]) / (years[i + 1] - years[i])
        return pops[i] + (pops[i + 1] - pops[i]) * frac
    end
end
interp_pop = DataFrame(Year=all_years, Population=pop_vec)
Plots.plot(interp_pop.Year, interp_pop.Population)


# Plots.plot(human_pop.Year, human_pop.Population)
# Plots.plot(sugar_cane.Year, sugar_cane.Area)
