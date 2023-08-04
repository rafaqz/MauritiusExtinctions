using CSV, DataFrames, IUCNRedList, JSON3, TerminalPager
using DimensionalData
using Statistics, StatsBase
using GLMakie
using TextAnalysis
using Plots
using StatsPlots
using DimensionalData.LookupArrays

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
timerange = 1440:20:2000
iucn_assesment.decade_last_seen = iucn_assesment.yearLastSeen_cleaned .÷ step(timerange) .* step(timerange)

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

cause_labels = [
    "Residential & commercial"
    "Agriculture & aquaculture"
    "Energy production & mining"
    "Transport corridors"
    "Biological resource use"
    "Human  disturbance"
    "Natural system modifications"
    "Invasive & diseases"
    "Pollution"
    "Geological events"
    "Climate & weather"
    "Other options"
]
function count_causes(codes)
    cs = map(1:12) do i
        count(==(i), last.(codes))
    end
    DimArray(cs, Dim{:cause}(DimensionalData.Categorical(cause_labels)))
end

flat_threats_with_assesment = leftjoin(flat_threats, iucn_assesment; on=:name => :scientificName)

group_flat_threats = map(groups) do group
    get_flat_threats(get_group_dict(group, iucn_threats_dict))
end
group_threats_with_iucn = map(group_flat_threats) do flat_threats
    leftjoin(flat_threats, iucn_assesment; on=:name => :scientificName)
end
threats_centuries = DimArray(Ti(1400:2000; sampling=Intervals(Start()))) do century
    map(group_threats_with_iucn) do group
        filter(group) do row
            ismissing(row.yearLastSeen_cleaned) ? false : (row.yearLastSeen_cleaned in century..century+99)
        end
    end
end
# allcodes = union(flat_threats_with_assesment.name .=> parse.(Int, first.(split.(flat_threats.code, '.'))))
group_codes = map(group_flat_threats) do flat_threats
    union(flat_threats.name .=> parse.(Int, first.(split.(flat_threats.code, '.'))))
end
group_causes = map(count_causes, group_codes)
cause_matrix = rebuild(cat(group_causes...; dims=Dim{:group}(collect(keys(group_causes)))); name=:count)
century_group_codes = map(threats_centuries) do century_groups
    map(century_groups) do flat_threats
        union(flat_threats.name .=> parse.(Int, first.(split.(flat_threats.code, '.'))))
    end
end
century_group_causes = map(century_group_codes) do group_codes
    map(count_causes, group_codes)
end
century_cause_matrix = map(century_group_causes) do group_causes
    rebuild(cat(group_causes...; dims=Dim{:group}(collect(keys(group_causes)))); name=:count)
end
# all_causes = count_causes(allcodes)
C = 1700
for C in 1400:100:2000 
cause_matrix = century_cause_matrix[Ti=At(C)]
# Bar plot of all causes by group
ps = map(1:12) do i
    bar(cause_matrix[cause=i]; ylims=(0, 120))
end
display(Plots.plot(ps...; size=(1900, 1900)))
sleep(1)
end

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

function get_threat_timeline(df)
    threat_timeline = map(x -> Int[], Ti(timerange))
    for (d, tc) in zip(df.decade_last_seen, df.threat_codes)
        (ismissing(tc) || ismissing(d)) && continue 
        a = threat_timeline[At(d)]
        append!(a, tc)
    end
    threat_timeline_counts = map(threat_timeline) do xs
        a = map(Dim{:cause}(DimensionalData.Categorical(1:12))) do i
            count(==(i), xs)
        end
        set(a, :cause => cause_labels)
    end
end

as_matrix(timeline) = cat(timeline...; dims=dims(timeline, Ti))
function get_threat_matrix(df)
    timeline = get_threat_timeline(df)
    threat_matrix = as_matrix(timeline)
    return threat_matrix
end
tps = map(collect(pairs(groups))) do (k, group)
    Plots.plot(get_threat_matrix(group); title="Causes: $k")
end;
Plots.plot(tps...; size=(1400, 1800), layout=(4, 2))
group_threat_timelines = map(get_threat_timeline, groups)
group_extinction_counts = map(by_decade, groups)
cause_freqs = map(group_threat_timelines, group_extinction_counts) do tt, ec
    as_matrix(tt ./ ec)
end
rtps = map(collect(pairs(cause_freqs))) do (k, m)
    Plots.plot(m; title="Cause frequency: $k", legend=:topleft, markerstrokewidth=0)
end;
Plots.plot(rtps...; size=(1400, 1800), layout=(4, 2))

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

by_decade(df) = by_decade(df) do d
    n = nrow(d)
    n == 0 ? missing : n
end
function by_decade(f, df)
    DimArray(Ti(timerange)) do decade
        sub = subset(df, :decade_last_seen => ByRow(==(decade)); skipmissing=true)
        f(sub)
    end
end
function uncertainty_decade(df)
    by_decade(df) do dc
        if nrow(dc) == 0
            missing
        else
            mean(dc.uncertainty_bool)
        end
    end
end

plotgroups = (:bird, :mammal, :reptile, :molusc, :fish, :plant)
all_extinction_decades = by_decade(iucn_assesment)
p1_kw = (palette=:Dark2_5, legend=:topleft, opacity=0.5, markersize=4, markerstrokewidth=0)
p1 = Plots.scatter(all_extinction_decades; label="all", color=:black, p1_kw...)
group_extinction_decades = map(by_decade, groups[plotgroups])
for (k, dc) in pairs(group_extinction_decades)
    Plots.scatter!(p1, dc; title="Recorded Extinctions", label=string(k), p1_kw...)
end
uncertainty_decades = map(uncertainty_decade, groups)
all_uncertainty_decades = uncertainty_decade(iucn_assesment)
p2_kw = (legend=:bottomcenter, lw=2)
p2 = Plots.plot(all_uncertainty_decades; label="all", color=:black, p2_kw...)
group_uncertainty_decades = map(uncertainty_decade, groups[plotgroups])
for (k, dc) in pairs(group_uncertainty_decades)
    Plots.plot!(p2, dc; label=string(k), title="Extinction Cause Uncertainty", p2_kw...)
end
Plots.plot(p1, p2; layout=(2, 1), size=(700, 1200))

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
     plot_kw=(; color=:lightblue, label="All", 
        title="Extinction Cause Uncertainty", 
        xguide="Cause uncertainty"
     ),
)
_bar(mean_text_uncertainty, groups; 
     mod=skipmissing, plot_kw=(; color=:blue, label="Threats\nassessed"), p=bar!,
)

schema_map = (
    hunting = ["8.1.*"],
    habitat = [""],
    invasives = [""],
    others = [""],
)
