using CSV, DataFrames, IUCNRedList, JSON3, TerminalPager
using DimensionalData
using Statistics, StatsBase
using GLMakie
using TextAnalysis
using Plots
using StatsPlots

# IUCNRedList.set_token("486f559a61285ba396234fc186897b94eda1bd15aa4216a8e1d9f4a8cf40d4c7")
# spec = species_by_category("EX") 

# threats = Dict{String,Any}[]
# foreach(iucn_assesment.scientificName) do name
#     println("\nFetching $(name)...")
#     threat = IUCNRedList.threats(String(name))
#     push!(threats, threat)
#     if !isempty(threat["result"])
#         println("Has result:")
#         display(threat["result"])
#     end
#     JSON3.write(iucn_threats_json_path, threats)
# end
# filter(threats) do threat
#     !isempty(threat["result"])
# end

iucn_threats_dict_json_path = "/home/raf/PhD/IUCN/iucn_threats_dict.json"
# iucn_threats_dict = Dict{String,Any}()
# for threat in threats
    # iucn_threats_dict[threat["name"]] = threat["result"]
# end
# JSON3.write(iucn_threats_dict_json_path, iucn_threats_dict)

iucn_assesment_path = "/home/raf/PhD/IUCN/iucn_assessments_gbif.csv"
iucn_assesment = CSV.read(iucn_assesment_path, DataFrame)
sort!(iucn_assesment, :scientificName)
pairs(names(iucn_assesment))
# iucn_threats_json_path = "/home/raf/PhD/IUCN/iucn_threats.json"
# run(`libreoffice $iucn_assesment_path`)
iucn_assesment.uncertainty_bool = map(iucn_assesment.uncertainty) do u
    ismissing(u) || u == "Yes"
end
iucn_assesment.decade_last_seen = iucn_assesment.yearLastSeen_cleaned .÷ 20 .* 20

# Remove html
cleaned_threats = map(iucn_assesment.threats) do threat
    if ismissing(threat) || length(threat) == 0
        missing
    else
        doc = StringDocument(threat)
        prepare!(doc, strip_html_tags)
        TextAnalysis.text(doc)
    end
end
iucn_assesment.cleaned_threats .= cleaned_threats
iucn_assesment.threats_length = map(iucn_assesment.cleaned_threats) do ct
    ismissing(ct) ? missing : length(ct)
end
#
iucn_threats_dict = JSON3.read(iucn_threats_dict_json_path, Dict{String,Any})
severity = map(collect(pairs(iucn_threats_dict))) do (k, v)
    k => map(v) do threat
        if haskey(threat, "severity")
            threat["severity"]
        else
            nothing
        end
    end
end |> Dict

severity_uncertainty = map(collect(pairs(severity))) do (k, v)
    length(v) == 0 && return k => 1.0
    k => mean(v) do x
        isnothing(x) || x == "Unknown"
    end
end |> sort
sort!(severity_uncertainty)
insertcols!(iucn_assesment, 14, :severity_uncertainty => last.(severity_uncertainty))
# CSV.write(iucn_assesment_path, iucn_assesment)
# Mean ceertainty
mean(last.(severity_uncertainty))
# How many have any certainty at all
c = count(<(1), last.(severity_uncertainty))
c / length(severity_uncertainty)

# Thread classifications
function get_flat_threats(threats_dict)
    threatkeys = nothing
    for (k, v) in threats_dict
        if length(v) > 0
            threatkeys = Tuple(Symbol.(keys(v[1]))) 
        end
    end
    isnothing(threatkeys) && return missing
    allkeys = (:name, threatkeys...)
    map(collect(pairs(threats_dict))) do (k, v)
        map(v) do threat
            as_missings = map(v -> isnothing(v) ? missing : v, values(threat))
            NamedTuple{allkeys}((k, ntuple(i -> as_missings[i], 7)...))
        end
    end |> Iterators.flatten |> collect |> DataFrame
end
flat_threats = get_flat_threats(iucn_threats_dict)

function get_group_dict(group, threats)
    group_threats = Dict{String,Any}()
    for name in group.scientificName
        group_threats[name] = get(threats, name, missing)
    end
    return group_threats
end
group_flat_threats = map(groups) do group
    get_flat_threats(get_group_dict(group, iucn_threats_dict))
end
group_flat_threats.amphibian
allcodes = union(flat_threats.name .=> parse.(Int, first.(split.(flat_threats.code, '.'))))
group_codes = map(group_flat_threats) do flat_threats
    union(flat_threats.name .=> parse.(Int, first.(split.(flat_threats.code, '.'))))
end

cause_labels = [
    "Residential & commercial development"
    "Agriculture & aquaculture"
    "Energy production & mining"
    "Transportation & service corridors"
    "Biological resource use"
    "Human intrusions & disturbance"
    "Natural system modifications"
    "Invasive & other problematic species, genes & diseases"
    "Pollution"
    "Geological events"
    "Climate change & severe weather"
    "Other options"
]
function count_causes(codes)
    cs = map(1:12) do i
        count(==(i), last.(codes))
    end
    DimArray(cs, Dim{:cause}(DimensionalData.Categorical(cause_labels)))
end
all_causes = count_causes(allcodes)
group_causes = map(count_causes, group_codes)
cause_matrix = rebuild(cat(group_causes...; dims=Dim{:group}(collect(keys(group_causes)))); name=:count)
# ctg = repeat(DimensionalData.lookup(cause_matrix, :group), size(cause_matrix, :cause))
# nam = repeat(DimensionalData.lookup(cause_matrix, :cause), size(cause_matrix, :group))
# cause_matrix[cause=8:8]
# groupedbar(nam, cause_matrix, group = ctg, xlabel = "Groups", ylabel = "Scores",
#     title="", bar_width=0.67, lw=0, framestyle=:box,
#     xticks= DimensionalData.lookup(cause_matrix, :cause)
# )

# Bar plot of all causes by group
ps = map(1:12) do i
    bar(cause_matrix[cause=i]; ylims=(0, 120))
end
Plots.plot(ps...; size=(1900, 1900))

grouped_codes = Dict{String,Any}()
for (k, v) in allcodes
    if haskey(grouped_codes, k)
        push!(grouped_codes[k], v)
    else
        grouped_codes[k] = [v]
    end
end
species_threat_codes = map(iucn_assesment.scientificName) do name
    get(grouped_codes, name, missing)
end
iucn_assesment.threat_codes = species_threat_codes
bru = filter(:threat_codes => c -> ismissing(c) ? false : 5 in c, iucn_assesment) 
sort!(iucn_assesment, :yearLastSeen_cleaned)
threat_timeline = map(x -> Int[], Ti(1440:20:2000))
for (d, tc) in zip(iucn_assesment.decade_last_seen, iucn_assesment.threat_codes)
    (ismissing(tc) || ismissing(d)) && continue 
    a = threat_timeline[At(d)]
    @show a tc
    append!(a, tc)
end
threat_timeline
threat_timeline_counts = map(threat_timeline) do xs
    a = map(Dim{:cause}(DimensionalData.Categorical(1:12))) do i
        count(==(i), xs)
    end
    set(a, :cause => cause_labels)
end
threat_matrix = cat(threat_timeline_counts...; dims=dims(threat_timeline_counts, Ti))
Plots.plot(threat_matrix)

flat_threats_path = "/home/raf/PhD/IUCN/flat_threats.csv"
# CSV.write(flat_threats_path, flat_threats);
# run(`libreoffice $flat_threats_path`)

class_filters = (
    arthropod = :phylumName => ByRow(==("ARTHROPODA")),
    amphibian = :className => ByRow(==("AMPHIBIA")),
    bird = :className => ByRow(==("AVES")),
    fish = :className => ByRow(x -> x in ("ACTINOPTERYGII", "SARCOPTERYGII", "EUSELACHII", "HOLOCEPHALI")),
    mammal = :className => ByRow(==("MAMMALIA")),
    molusc = :phylumName => ByRow(==("MOLLUSCA")),
    plant = :kingdomName => ByRow(==("PLANTAE")),
    reptile = :className => ByRow(==("REPTILIA")),
)
groups = map(class_filters) do filter
    subset(iucn_assesment, filter) 
end
group_sizes = map(nrow, groups) |> pairs

by_decade(df) = by_decade(nrow, df)
function by_decade(f, df)
    map(Ti(1500:20:2000)) do decade
        sub = subset(df, :decade_last_seen => ByRow(==(decade)); skipmissing=true)
        f(sub)
    end
end

plotgroups = (:bird, :mammal, :reptile, :molusc, :fish, :plant)
all_extinction_decades = by_decade(iucn_assesment)
p1 = Plots.scatter(all_extinction_decades; palette=:Dark2_5, markersize=3, label="all", color=:black)
group_extinction_decades = map(by_decade, groups[plotgroups])
for (k, dc) in pairs(group_extinction_decades)
    Plots.scatter!(p1, dc; 
        title="Extinctions", 
        label=string(k), markersize=3, opacity=0.5
    )
end
p1

function uncertainty_decade(df)
    by_decade(df) do dc
        if nrow(dc) == 0
            missing
        else
            mean(dc.uncertainty_bool)
        end
    end
end
uncertainty_decades = map(uncertainty_decade, groups)
all_uncertainty_decades = uncertainty_decade(iucn_assesment)
p2 = Plots.plot(all_uncertainty_decades; 
    palette=:Dark2_5, markersize=3, label="all", color=:black
)
group_uncertainty_decades = map(uncertainty_decade, groups[plotgroups])
for (k, dc) in pairs(group_uncertainty_decades)
    Plots.plot!(p2, dc; 
        title="Extinction Cause Uncertainty", 
        label=string(k), markersize=3, opacity=0.5, legend=:bottomleft
    )
end
p2

Plots.plot(p1, p2; layout=(2, 1), size=(700, 1200))

p = Plots.scatter(first.(all_extinction_decades), last.(all_extinction_decades); opacity=0.5)

mean_severity_uncertainty(f, x) = mean(f(x.severity_uncertainty))
mean_citations(f, x) = mean(f(x.citation .== Ref("Yes")))
mean_threat_length(f, x) = mean(f(x.threats_length)) / 1000
mean_uncertainty_and_citations(f, x) = mean(f(x.severity_uncertainty .* x.citation .== Ref("Yes")))
mean_text_uncertainty(f, x) = mean(f(x.uncertainty .== Ref("Yes")))

function _bar(f, table; mod=identity, plot_kw=(;), p=Plots.bar) 
    xs = map(t -> f(mod, t), table) |> pairs |> collect
    p(map(last, xs);
        orientation=:h, yticks=(1:8, map(string ∘ first, xs)), xlims=(0, 1),
        plot_kw...,
    )
end
_bar(mean_text_uncertainty, groups; 
     mod=xs->replace(xs, missing => 1.0),
     plot_kw=(; color=:lightblue, label="All", title="Extinction Cause Uncertainty"),
)
_bar(mean_text_uncertainty, groups; 
     mod=skipmissing, plot_kw=(; color=:blue, label="Assessed"), p=bar!,
)

schema_map = (
    hunting = ["8.1"],
    habitat = [""],
    invasives = [""],
    others = [""],
)
