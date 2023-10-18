function filter_population(table)
    periods = names(table)[2:end-1]
    times = map(periods) do p
        s = split(p, '-')
        Base.parse(Int, s[1]), Base.parse(Int, s[2])
    end
    boundsmatrix = reinterpret(reshape, Int, times)
    timedim = Ti(Sampled(map(first, times);
        order=ForwardOrdered(),
        span=Explicit(boundsmatrix),
        sampling=Intervals(Start())),
    )
    populations = Dict{Symbol,Any}()
    for i in 1:size(table, 1)
        popvals = Vector{Union{Missing,Int}}(undef, length(periods))
        popvals .= missing
        started = false
        rawdata = table[i, 2:end-1]
        common_name = table[i, 1]
        local lastval = missing
        for j in eachindex(popvals)
            !started && ismissing(rawdata[j]) && continue
            rawval = rawdata[j]
            rawval = rawval isa String ? rawval[1:1] : rawval
            if haskey(data_key, rawval)
                category = data_key[rawval]
            else
                @warn "unidentified value in table: $rawval, for $common_name. `missing` used instead"
                continue
            end
            pop_estimate = if ismissing(category)
                lastval
            elseif category == "abundant" 
                3
            elseif category == "common" 
                2
            elseif category == "rare" 
                1
            elseif category == "uncertain" 
                0
                # ismissing(lastval) ? 1 : lastval
                # lastval
            elseif category == "extinction/absence"
                0
            elseif category == "observed"
                ismissing(lastval) ? 1 : lastval
            elseif category == "likely but unconfirmed"
                1
            elseif category == "several species unseparated"
                ismissing(lastval) ? 1 : lastval
            elseif category == "present no record"
                lastval
            elseif category == "not reported"
                lastval
            elseif category == "recorded introduction"
                1
            elseif category == "approximate introduction"
                1
            elseif category == "captive only"
                0
            elseif category == "move out of category" 
                lastval
                # @info "category switch in $(table[i, 1])"
                # :moveout
            elseif category == "move into category"
                lastval
                # @info "category switch in $(table[i, 1])"
                # :movein
            else
                @show "Unfound: $category"
                continue
            end
            lastval = pop_estimate
            if ismissing(popvals[j])
                popvals[j] = pop_estimate
            end
            started = true
        end
        for i in reverse(eachindex(popvals))
            if ismissing(rawdata[i])
                popvals[i] = missing
            end
        end
        populations[Symbol(common_name)] = DimArray(popvals, timedim)
    end
    return DataFrame(:period => periods, populations...)
end

reverse_data_key = OrderedDict(
    "abundant" => "a",
    "common" => "b",
    "rare" => "c",
    "uncertain" => "d",
    "extinction/absence" => "e",
    "observed" => "f",
    "likely but unconfirmed" => "g",
    "several species unseparated" => "h",
    "present no record" => "i",
    "not reported" => "j",
    "recorded introduction" => "k",
    "approximate introduction" => "L",
    "captive only" => "m",
    "move out of category"  =>  "N",
    "move into category" => "O",
    missing => missing,
)
data_key = OrderedDict(reverse(p) for p in reverse_data_key)

as_dataframe(xl, name::Symbol) = as_dataframe(xl, sheetnames[name])
function as_dataframe(xl, name::String)
    sheet = xl[name]
    df = DataFrame(XLSX.eachtablerow(sheet))
    df.Species .= strip.(isnumeric, df.Species)
    return df
end

function filter_observations(table, species_lookup)
    periods = names(table)[2:end-1]
    times = map(periods) do p
        s = split(p, '-')
        Base.parse(Int, s[1]), Base.parse(Int, s[2])
    end
    boundsmatrix = reinterpret(reshape, Int, times)
    timedim = Ti(Sampled(map(first, times);
        order=ForwardOrdered(),
        span=Explicit(boundsmatrix),
        sampling=Intervals(Start())),
    )
    species = Dict{Symbol,Any}()
    for i in 1:size(table, 1)
        common_name = table[i, 1]
        haskey(species_lookup, common_name) || continue
        observations = falses(length(periods))
        rawdata = table[i, 2:end-1]

        for j in eachindex(observations)
            rawval = rawdata[j]
            rawval = rawval isa String ? rawval[1:1] : rawval
            if haskey(data_key, rawval)
                category = data_key[rawval]
            else
                @warn "unidentified value in table: $rawval, for $common_name. `missing` used instead"
                continue
            end
            present = if ismissing(category)
                false
            elseif category == "abundant" 
                true
            elseif category == "common" 
                true
            elseif category == "rare" 
                true
            elseif category == "uncertain" 
                0
                # ismissing(lastval) ? 1 : lastval
                # lastval
            elseif category == "extinction/absence"
                false
            elseif category == "observed"
                true
            elseif category == "likely but unconfirmed"
                1
            elseif category == "several species unseparated"
                true
                # ismissing(lastval) ? 1 : lastval
            elseif category == "present no record"
                false
            elseif category == "not reported"
                false
            elseif category == "recorded introduction"
                true
            elseif category == "approximate introduction"
                true
            elseif category == "captive only"
                false
            elseif category == "move out of category" 
                false
            elseif category == "move into category"
                false
            else
                @show "Unfound: $category"
                continue
            end
            observations[j] = present
        end
        species[species_lookup[common_name]] = DimArray(observations, timedim)
    end
    return DataFrame(:period => periods, species...)
end
