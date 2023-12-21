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
rod_pop = CSV.File(joinpath(datadir, "Population/rodrigues.csv")) |> DataFrame

human_pop_timelines = map((mus=mus_pop, reu=reu_pop, rod=rod_pop)) do pop
    pop_vec = interpolate_years(pop, (:Year, :Population))
    DimArray(pop_vec, Ti(all_years); name=:Human_Population)
end
human_pop_timelines.rod[At(1810)]

sugar_vec = interpolate_years(sugar_cane, (:Year, :Area))
sugar_timeline = DimArray(sugar_vec, Ti(all_years); name=:Area_Sugar)


# using CairoMakie

# set_theme!(theme_ggplot2())
# Makie.scatter()
# Makie.scatter!(sugar_timeline .* 10000)
# fig = Figure(resolution = (400, 400));
# ax1 = fig[1:5, 1:4] = Axis(fig, xlabel = "year")
# ax2 = fig[1:5, 1:4] = Axis(fig)

# plot_kw = (; transparency=true, linewidth=2)
# plot1 = Makie.lines!(ax1, lookup(human_pop_timeline.mus, Ti), human_pop_timeline.mus; 
#     color=(:lightblue, 0.6), plot_kw...
# )
# plot2 = Makie.lines!(ax2, lookup(sugar_timeline, Ti), sugar_timeline; 
#     color=(:red, 0.6), plot_kw...
# )
# ax2.yaxisposition = :right
# ax2.yticklabelalign = (:left, :center)
# ax2.xticklabelsvisible = false
# ax2.xticklabelsvisible = false
# ax1.ygridvisible=false
# ax2.ygridvisible=false
# ax2.xlabelvisible = false
# linkxaxes!(ax1,ax2)
# Legend(fig[6, 1:4], [plot1, plot2], ["Human Population", "Sugarcane expansion"])
# fig
# save("lucc_drivers.png", fig)

# Plot
# Plots.plot(human_pop_timeline.mus)
# Plots.plot(sugar_timeline)
