using GLMakie, Unitful
using CairoMakie
using DataFrames, Dates
CairoMakie.activate!()
# GLMakie.activate!()

include("species_tables.jl")
include("tabular_data.jl")
include("plots/original_vegetation.jl")

simple_growth(N, k, r, t) = (N .* k) ./ (N .+ (k.- N) .* exp.(.-(r * t)))

human_habitation_periods = map(human_pop_timelines) do A
    A[findfirst(>(0), A):end]
end
labels = ["Mauritius", "Reunion", "Rodrigues"]

bird_energy_content = 10.0u"kJ/g" # roughly
areas = (mus=2040u"km^2", reu=2512u"km^2", rod=108u"km^2")
kcats = 1u"km^-2"
krats = 30u"ha^-1"

maximum(island_endemic_tables.mus.Mass)
maximum(island_endemic_tables.rod.Mass)

birds_by_humans = nhumans * 0.2u"d^-1" * 365u"d" / 2040u"km^2"

birds_by_cats = ncats * 0.02u"d^-1" * 365u"d" / 2040u"km^2"

max_cats = 1u"km^-2" * 2040u"km^2" # ?
max_rats = 30u"ha^-1" * 2040u"km^2" # ?
cat_introduced = 1688

human_energy_intake = 8700u"kJ/d"
cat_energy_intake = 2000u"kJ/d"
rat_energy_intake = 350u"kJ/d"
pig_energy_intake = 9000u"kJ/d" # Made up, get the real numbers
assimilation_efficiency = 0.84
human_bird_diet_fraction = 0.1
cat_bird_diet_fraction = 0.1
rat_bird_diet_fraction = 0.0001
pig_bird_diet_fraction = 0.01

introductions = map(island_keys) do k
    map((cat="cat", black_rat="black_rat", pig="pig")) do  x
        subset(introductions_df, :species => ByRow(==(x)), :island => ByRow(==(string(k)))).year[1]
    end
end

extinct = subset(island_endemic_tables.mus, :extinct => ByRow(!ismissing))
extinct = subset(island_endemic_tables.reu, :extinct => ByRow(!ismissing))
extinct = subset(island_endemic_tables.rod, :extinct => ByRow(!ismissing))

human_bird_consumption = human_bird_diet_fraction * human_energy_intake / (bird_energy_content * assimilation_efficiency)
cat_bird_consumption = cat_bird_diet_fraction * cat_energy_intake / (bird_energy_content * assimilation_efficiency)
rat_bird_consumption = rat_bird_diet_fraction * rat_energy_intake / (bird_energy_content * assimilation_efficiency)
pig_bird_consumption = pig_bird_diet_fraction * pig_energy_intake / (bird_energy_content * assimilation_efficiency)

eaten_by_humans = human_pop_timelines.mus[1600 .. 1720] .* human_bird_consumption
eaten_by_cats = max_cats * length((cat_introduced+10):1720) *  human_bird_consumption
eaten_by_rats = max_rats * length(1600:1720) * rat_bird_consumption

human_consumption_rate = uconvert(u"kg/yr*km^(-2)", mean(eaten_by_humans) / 2040u"km^2")
cat_consumption_rate = uconvert(u"kg/yr*km^(-2)", mean(eaten_by_cats) / 2040u"km^2")
rat_consumption_rate = uconvert(u"kg/yr*km^(-2)", mean(eaten_by_rats) / 2040u"km^2")

eaten_by_humans = map(human_habitation_periods) do human_pops
    uconvert.(u"g/d", human_pops .* human_bird_consumption)
end
t = 1u"yr"/12
plotstop = 2200
eaten_by_cats = map(introductions, areas) do (; cat), area
    r = eachrow(pred_df)[1].rmax * u"yr^-1"
    global P = 10
    cat_pops = DimArray(Ti(cat:1//12:plotstop)) do d
        P = simple_growth(P, 0.01u"ha^-1" * area, r, t)
    end
    uconvert.(u"g/d", cat_pops .* cat_bird_consumption)
end
eaten_by_rats = map(introductions, areas) do (; black_rat), area
    r = eachrow(pred_df)[3].rmax * u"yr^-1"
    global P = 10
    rat_pops = DimArray(Ti(black_rat:1//12:plotstop)) do d
        P = simple_growth(P, 30u"ha^-1" * area, r, t)
    end
    uconvert.(u"g/d", rat_pops .* rat_bird_consumption)
end
eaten_by_pigs = map(introductions, areas) do (; pig), area
    r = eachrow(pred_df)[5].rmax * u"yr^-1"
    global P = 10
    pig_pops = DimArray(Ti(pig:1//12:plotstop)) do d
        P = simple_growth(P, 0.02u"ha^-1" * area, r, t)
    end
    uconvert.(u"g/d", pig_pops .* pig_bird_consumption)
end

function plot_consumption!(ax, by_cats, by_rats, by_pigs, by_humans, i; fontsize=12)
    Makie.lines!(ax, lookup(by_cats, Ti), ustrip.(u"g/d", by_cats); color=:orange, label="Cats 10% bird diet")
    Makie.lines!(ax, lookup(by_cats, Ti), ustrip.(u"g/d", by_cats .* 5); color=(:orange, 0.5), label="Cats 50% bird diet")
    Makie.lines!(ax, lookup(by_rats, Ti), ustrip.(u"g/d", by_rats); color="#067dd1", label="Rats 0.01% bird diet")
    Makie.lines!(ax, lookup(by_rats, Ti), ustrip.(u"g/d", by_rats .* 5); color=("#067dd1", 0.5), label="Rats 0.02% bird diet")
    Makie.lines!(ax, lookup(by_humans, Ti), ustrip.(u"g/d", by_humans); color=:black, label="Human 10% bird diet")
    Makie.lines!(ax, lookup(by_humans, Ti), ustrip.(u"g/d", by_humans ./ 500); color=(:black, 0.5), label="Human 0.02% bird diet")
    # i == 1 && axislegend(ax; position=:lt)
    if length(by_pigs) > 0 && by_pigs[Near(1650)] > 0u"g/d"
        Makie.lines!(ax, lookup(by_pigs, Ti), ustrip.(u"g/d", by_pigs); color=:green, label="Pig 1% bird diet")
        Makie.text!(ax, 1660, ustrip(u"g/d", by_pigs[Near(1660)]); text="Pig 1% bird diet", fontsize)
    end
    Makie.text!(ax, 1600, ustrip(u"g/d", by_rats[Near(1600)]); text="Black rat 0.01% bird diet", fontsize)
    Makie.text!(ax, 1600, ustrip(u"g/d", by_rats[Near(1600)] * 5); text="Black rat 0.05% bird diet", fontsize)
    Makie.text!(ax, 1750, ustrip(u"g/d", by_cats[Near(1750)]); text="Cat 10% bird diet", fontsize)
    Makie.text!(ax, 1750, ustrip(u"g/d", by_cats[Near(1750)] * 5); text="Cat 50% bird diet", fontsize)
    Makie.text!(ax, 1741, ustrip(u"g/d", by_humans[Near(1739)]); text="Human 10% bird diet", fontsize)
    Makie.text!(ax, 1980, ustrip(u"g/d", by_humans[Near(1980)] / 500); align=(:right, :center), text="Human 0.05% bird diet", fontsize)
end

function extinction_plot!(ax, table; 
    label, 
    textcolumn=:LostLand_name, 
    fontsize=10,
    color=:black,
    text=true,
) 
    extinct = subset(table, :extinct => ByRow(!ismissing))
    extinct.extinct1 = map(extinct.extinct) do r
        if last(r) == 2100 
            first(r) - 2.5, first(r) + 2.5
        else
            first(r), last(r)
        end
    end
    # period_len = 50
    # extinct.decade = round.(Int, (mean.(extinct.extinct1))) .÷ period_len .* period_len
    # period_mean_mass = map(df -> mean(df.Mass), collect(DataFrames.groupby(extinct, :decade; sort=true)))
    # periods = sort!(union(extinct.decade)) .+ (period_len / 2)
    # hunted = subset(extinct, :Hunting_preference => ByRow(>(0)))
    # p = Makie.scatter(hunted.extinct1, hunted.Mass; markersize=30, color=(:red, 0.2), axis=())
    segments = map(extinct.extinct1, extinct.Mass) do (x1, x2), y
        (Makie.Point(x1, y), Makie.Point(x2, y)) 
    end
    Makie.linesegments!(ax, segments; color, label, linewidth=1.5)
    # Makie.lines!(ax, periods, period_mean_mass; color=(color, 0.5), linewidth=1, label="Mean 50yr mass")
    Makie.hlines!(ax, median(extinct.Mass); color=(color, 0.5), linewidth=1, linestyle=:dash, label="Median mass")
    if text
        Makie.text!(ax, mean.(extinct.extinct1), extinct.Mass .* 1.15; 
            text=map(s -> ismissing(s) ? "" : s, getproperty(extinct, textcolumn)), # .* '\n' .* extinct.Species)
            align=(:center, :baseline),
            fontsize,
        )
    end
end

###############################################################################################
# Whole page plot
fig = Figure(; size=(800, 800), title="Extinction times by Mass");
weights = [1, 10, 100, 1000, 10000]
axs1 = map(1:3) do i
    ax = Axis(fig[i * 3 - 2, 2]; 
        xlabel="Year", 
        ylabel="Extinction\nby Mass (g)",
        yscale=log10,
        xticks=1600:100:2000
    )
    i == 3 || hidexdecorations!(ax; grid=false)
    Makie.ylims!(ax, (10, 280000))
    Makie.hidespines!(ax)
    ax
end
axs2 = map(1:2) do i
    ax = Axis(fig[i * 3 - 1, 2]; 
        xlabel="Year", 
        ylabel="Habitat class",
        xticks=1600:100:2000
    )
    i == 3 || hidexdecorations!(ax; grid=false)
    Makie.ylims!(ax, (0.00000000000001, 1))
    Makie.hidespines!(ax)
    ax
end
axs2b = map(1:3) do i
    ax = Axis(fig[i * 3 - 1, 2]; 
        ylabel="Human population",
        yaxisposition = :right,
        xticks=1600:100:2000
    )
    i == 3 || hidexdecorations!(ax; grid=false)
    hidespines!(ax)
    ax
end
axs3 = map(1:3) do i
    ax = Axis(fig[i * 3, 2];
        xlabel="Year", 
        ylabel="Bird mass\nconsumed (kg)",
    )
    Makie.ylims!(ax, (10, 3e5))
    i == 3 || hidexdecorations!(ax; grid=false)
    hidedecorations!(ax; label=false)
    hidespines!(ax)
    ax
end
cmap = :tableau_20
habitat_colors = (
    mus=map(x -> getproperty(ColorSchemes, cmap)[(x - 1) / 9 ], 1:nhabitats.mus) |> reverse,
    reu=map(x -> getproperty(ColorSchemes, cmap)[(x - 1) / (nhabitats.reu - 1) ], 1:nhabitats.reu),
)
foreach(plot_aggregate!, axs2, uncleared, habitat_colors)
linkxaxes!(axs1..., axs2..., axs2b..., axs3...)
Makie.xlims!(axs1[1], (1600, 2020))
foreach(axs1, island_endemic_tables, labels) do ax, tbl, label
    extinction_plot!(ax, tbl; label="Extinction date range", color=:black)
    axislegend(ax; framevisible=false, labelsize=11)
end
foreach(axs2b, human_habitation_periods) do ax, human_pops
    Makie.lines!(ax, lookup(human_pops, 1), human_pops; color=(:black, 0.5))
    Makie.vlines!(ax, first(lookup(human_pops, 1)); color=(:black, 0.5))
end
foreach(1:3, labels) do i, ylabel
    ax = Axis(fig[i * 3 - 1:i * 3, 1]; ylabel, ylabelsize=20)
    hidedecorations!(ax; label=false)
    hidespines!(ax)
end
foreach(axs3, eaten_by_cats, eaten_by_rats, eaten_by_pigs, eaten_by_humans, 1:3) do args...
    plot_consumption!(args...)
end
colsize!(fig.layout, 1, 10)
Makie.xlims!(axs1[3], (1600, 2025))
Makie.xlims!(axs2[2], (1600, 2025))
Makie.xlims!(axs2b[3], (1600, 2025))
display(fig)
save("images/extinction_year_mass.png", fig)

###############################################################################################
# Separate slide plots
figs = map(title -> Figure(; size=(1200, 600), title), labels) 
axs1 = map(figs) do fig
    ax = Axis(fig[1, 1]; 
        xlabel="Year", 
        ylabel="Extinction\nby Mass (g)",
        yscale=log10,
        xticks=1600:100:2000
    )
    Makie.ylims!(ax, (10, 280000))
    Makie.hidespines!(ax)
    ax
end
axs2 = map(figs) do fig
    ax = Axis(fig[3, 1]; 
        xlabel="Year", 
        ylabel="Habitat class",
        xticks=1600:100:2000
    )
    Makie.ylims!(ax, (0.00000000000001, 1))
    Makie.hidespines!(ax)
    ax
end
axs2b = map(figs) do fig
    ax = Axis(fig[3, 1]; 
        ylabel="Human population",
        yaxisposition = :right,
        xticks=1600:100:2000
    )
    hidespines!(ax)
    ax
end
axs3 = map(figs, [4e5, 4e5, 2e4]) do fig, ymax
    ax = Axis(fig[2, 1];
        xlabel="Year", 
        ylabel="Mass\nconsumed (kg)",
    )
    Makie.ylims!(ax, (10, ymax))
    hidexdecorations!(ax; grid=false, label=false)
    hidespines!(ax)
    ax
end
foreach(plot_aggregate!, axs2, uncleared, habitat_colors)
linkxaxes!.(axs1, axs2, axs2b, axs3)
Makie.xlims!(axs1[1], (1600, 2020))
foreach(axs1, island_endemic_tables, labels) do ax, tbl, label
    extinction_plot!(ax, tbl; label="Extinction date range", color=:red, text=false)
    axislegend(ax; framevisible=false, labelsize=11)
end
foreach(axs2b, human_habitation_periods) do ax, human_pops
    Makie.lines!(ax, lookup(human_pops, 1), human_pops; color=(:black, 0.5))
    Makie.vlines!(ax, first(lookup(human_pops, 1)); color=(:black, 0.5))
end
foreach(axs3, eaten_by_cats, eaten_by_rats, eaten_by_pigs, eaten_by_humans, 1:3) do args...
    plot_consumption!(args...)
end
Makie.xlims!.(axs1, ((1600, 2025),))
Makie.xlims!.(axs2, ((1600, 2025),))
Makie.xlims!.(axs3, ((1600, 2025),))
display(figs[1])

display(figs[2])
display(figs[3])
map(island_keys, figs) do key, fig
    save("images/extinction_year_mass_$key.png", fig)
end


################################################################

fig = Figure(; size=(800, 800));
axs1 = map(1:3) do i
    ax = Axis(fig[i * 2 - 1, 1];
        ylabel="Extinct Species\nMass (g)",
        yscale=log10,
        xticks=1600:50:2000
    )
    hidexdecorations!(ax; grid=false)
    hidespines!(ax)
    Makie.xlims!(ax, (1590, 1750))
    ax
end
axs2 = map(1:3) do i
    ax = Axis(fig[i * 2, 1];
        ylabel="Mass\nconsumed (kg)",
        xticks=1600:50:2000
    )
    i == 3 || hidexdecorations!(ax; grid=false)
    hidespines!(ax)
    Makie.xlims!(ax, (1590, plotstop))
    Makie.ylims!(ax, (10, 280000))
    ax
end
foreach(axs1, island_endemic_tables, labels) do ax, tbl, label
    extinction_plot!(ax, tbl; label)
end
foreach(axs2, eaten_by_cats, eaten_by_rats, eaten_by_pigs, eaten_by_humans, 1:3) do args...
    plot_consumption!(args...)
end
# Makie.ylims!(axs1, (10, 1000))
linkxaxes!(axs1..., axs2...)
display(fig)

# GLMakie.activate!()
# CairoMakie.activate!()
fig = Figure(; size=(800, 1600));
weights = [1, 10, 100, 1000, 10000]
axs1 = map(1:3) do i
    ax = Axis(fig[i*2 - 1, 1]; 
        ylabel="Mass (g)",
        yscale=log10,
    )
    hidexdecorations!(ax; grid=false)
    ax
end
axs2 = map(1:3) do i
    ax = Axis(fig[i * 2, 1];
        xlabel="Year", 
        ylabel="Mass\nconsumed (kg)",
    )
    i == 3 || hidexdecorations!(ax; grid=false)
    ax
end
linkxaxes!(axs1..., axs2...)
Makie.xlims!(axs1[1], (1590, 1750))
Makie.ylims!(axs1[1], (10, 260000))
foreach(axs1, island_endemic_tables) do ax, tbl
    extinction_plot!(ax, tbl; label="Mauritius")
end
foreach(axs2, human_pop_timelines) do ax, human_pops
    Makie.lines!(ax, lookup(human_pops, 1), human_pops)
end
display(fig)

# Makie.axislegend(ax2; position=:lt)
