include("species_tables.jl")

mus_extinct = subset(island_endemic_tables.mus, :extinct => ByRow(!ismissing))
mus_extinct.logMass = log.(mus_extinct.Mass)
mus_extinct.extinct1 = first.(mus_extinct.extinct)
mus_hunted = subset(mus_extinct, :Hunting_preference => ByRow(>(0)))
p = Makie.scatter(mus_hunted.extinct1, mus_hunted.logMass; markersize=30, color=(:red, 0.4), axis=(xlabel="Year", ylabel="Log mass"))
Makie.scatter!(mus_extinct.extinct1, mus_extinct.logMass)
Makie.text!(mus_extinct.extinct1, mus_extinct.logMass; text=mus_extinct.Common_name .* '\n' .* mus_extinct.Species)
Makie.hlines!(mean(mus_extinct.logMass))
sort!(mus_extinct, :extinct1)
select(mus_extinct, :Common_name, :extinct1, :Mass)

reu_extinct = subset(island_endemic_tables.reu, :extinct => ByRow(!ismissing))
reu_extinct.logMass = log.(reu_extinct.Mass)
reu_extinct.extinct1 = first.(reu_extinct.extinct)
reu_hunted = subset(reu_extinct, :Hunting_preference => ByRow(>(0)))
Makie.scatter!(reu_hunted.extinct1, reu_hunted.logMass; markersize=30, color=(:red, 0.4))
Makie.scatter!(reu_extinct.extinct1, reu_extinct.logMass)
Makie.text!(reu_extinct.extinct1, reu_extinct.logMass; text=reu_extinct.Common_name .* '\n' .* reu_extinct.Species)
Makie.hlines!(mean(reu_extinct.logMass))
sort!(reu_extinct, :extinct1)
select(reu_extinct, :Common_name, :extinct1, :Mass)

rod_extinct = subset(island_endemic_tables.rod, :extinct => ByRow(!ismissing))
rod_extinct.Mass
rod_extinct.logMass = log.(rod_extinct.Mass)
rod_extinct.extinct1 = first.(rod_extinct.extinct)
rod_hunted = subset(rod_extinct, :Hunting_preference => ByRow(>(0)))
Makie.scatter!(rod_hunted.extinct1, rod_hunted.logMass; markersize=30, color=(:red, 0.4))
Makie.scatter!(rod_extinct.extinct1, rod_extinct.logMass)
Makie.text!(rod_extinct.extinct1, rod_extinct.logMass; text=rod_extinct.Common_name .* '\n' .* rod_extinct.Species)
Makie.hlines!(mean(rod_extinct.logMass))
sort!(rod_extinct, :extinct1)
select(rod_extinct, :Common_name, :extinct1, :Mass)



bird_energy_content = 10.0u"kJ/g" # roughly
maximum(island_endemic_tables.mus.Mass)
maximum(island_endemic_tables.rod.Mass)

birds_by_humans = nhumans * 0.2u"d^-1" * 365u"d" / 2040u"km^2"

birds_by_cats = ncats * 0.02u"d^-1" * 365u"d" / 2040u"km^2"

max_cats = 1u"km^-2" * 2040u"km^2" # ?
max_rats = 50u"ha^-1" * 2040u"km^2" # ?
cat_introduced = 1688

human_energy_intake = 8700u"kJ/d"
cat_energy_intake = 2131u"kJ/d"
rat_energy_intake = uconvert(u"kJ/d", 309u"kcal/d")
assimilation_efficiency = 0.84
human_bird_diet_fraction = 0.1
cat_bird_diet_fraction = 0.1
rat_bird_diet_fraction = 0.0001

human_bird_consumption = human_bird_diet_fraction * human_energy_intake / (bird_energy_content * assimilation_efficiency)
cat_bird_consumption = cat_bird_diet_fraction * cat_energy_intake / (bird_energy_content * assimilation_efficiency)
rat_bird_consumption = rat_bird_diet_fraction * rat_energy_intake / (bird_energy_content * assimilation_efficiency)
eaten_by_humans = human_pop_timelines.mus[1600 .. 1720] .* human_bird_consumption
eaten_by_cats = max_cats * length((cat_introduced+10):1720) *  human_bird_consumption
eaten_by_rats = max_rats * length(1600:1720) * rat_bird_consumption
human_consumption_rate = uconvert(u"kg/yr*km^(-2)", mean(eaten_by_humans) / 2040u"km^2")
cat_consumption_rate = uconvert(u"kg/yr*km^(-2)", mean(eaten_by_cats) / 2040u"km^2")
rat_consumption_rate = uconvert(u"kg/yr*km^(-2)", mean(eaten_by_rats) / 2040u"km^2")

ax2 = Axis(p.figure[2, 1])
Makie.scatter!(ax2, lookup(eaten, Ti), parent(ustrip.(eaten)))
linkxaxes!(p.axis, ax2)

