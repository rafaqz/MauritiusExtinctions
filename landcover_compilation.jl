# Rodrigues 50% cleared 1874 , Lost land page 151

# From woods and forests in Mauritius, p 40:40
# 1880: 70,000 acres out of 300,000 remain
# 35,000 were native (but "dilatipated and ruined?")
# Also mentions that invasives replace natives
include("travel_cost.jl")
include("tabular_data.jl")
include("map_file_functions.jl")
include("landcover_shared.jl")

include("map_file_list.jl")
compiled = compile_timeline(define_map_files(), masks, keys(lc_categories))

combined = map(compiled) do island
    Rasters.combine(namedvector_raster.(island.timeline))
end
merged = map(combined) do A
    cross_validate_timeline(A, transitions)
end
added = map(combined, merged) do raw, final
    map(raw, final) do rs, fs
    # We minimise forced values in the source throught its possible 
    # transitions, they may resolve to one or two distinct timelines. 
        if any(rs) # Ignore completely missing
            map(rs, fs) do r, f
                !r & f
            end
        else
            rs
        end
    end
end
filled = map(combined, merged) do raw, final
    map(raw, final) do rs, fs
        if any(rs)
            zero(rs)
        else
            fs
        end
    end
end
removed = map(combined, merged) do raw, final
    map(raw, final) do rs, fs
        map(rs, fs) do r, f
            r & !f
        end
    end
end
uncertain = map(merged) do final
    map(final) do fs
        if count(fs) > 1
            fs
        else
            zero(fs)
        end
    end
end
mixed_stats = (; combined, filled, uncertain, added, removed, merged)
landcover_statistics = map(mixed_stats...) do args...
    NamedTuple{keys(mixed_stats)}(args)
end
