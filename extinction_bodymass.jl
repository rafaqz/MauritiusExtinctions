include("species_tables.jl")

mus_extinct = subset(island_endemic_tables.mus, :extinct => ByRow(!ismissing))
mus_extinct.logMass = log.(mus_extinct.Mass)
mus_extinct.extinct1 = first.(mus_extinct.extinct)
mus_hunted = subset(mus_extinct, :Hunting_preference => ByRow(>(0)))
Makie.scatter(mus_hunted.extinct1, mus_hunted.logMass; markersize=30, color=(:red, 0.4), axis=(xlabel="Year", ylabel="Log mass"))
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

nhumans = 100 
birds_by_humans = nhumans * 0.2u"d^-1" * 365u"d" / 2040u"km^2"

ncats = 2040
birds_by_cats = ncats * 0.02u"d^-1" * 365u"d" / 2040u"km^2"

human_bird_diet_fraction = 0.1
human_energy_needs = 8700u"kJ/d"
human_bird_consumption = human_bird_diet_fraction * human_energy_needs / bird_energy_content
eaten = human_pop_timelines.mus[1600 .. 1730] .* human_bird_consumption
upreferred(sum(eaten .* u"d"))
Makie.scatter(1600:1730, ustrip.(eaten))

