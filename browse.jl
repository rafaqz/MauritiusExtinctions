
# Google search all species, opening the next when the browser is closed.
for row in eachrow(subset(species, :Origin => ByRow(==("Endemic"))))
    ismissing(row["Common name"]) || ismissing(row["Species"]) && continue
    placenamed = false
    for place in ("RéuReunion", "Mauritius", "Rodrigues", "Mascarene")
        placenamed = occursin(place, row["Common name"]) 
        @show row["Common name"] place placenamed
        placenamed && break
    end
    placenamed && continue
    search = replace(row["Species"], " " => "+")
    run(`chromium https\://www.google.com/search\?q=$(search)`)
    # run(`chromium https\://www.google.com/search\?tbm=isch\&q=$(search)`)
end

for row in Tables.rows(species)
    sp = row.Species
    row.Origin == "Alien" || continue
    ismissing(sp) && continue
    search = replace(sp, " " => "+")
    run(`chromium https\://www.google.com/search\?q=$(search)+pet`)
    # run(`chromium https\://www.google.com/search\?tbm=isch\&q=$(search)`)
end

# Check in Hume
humepath = "/home/raf/zotero_library/Journal Article/undefined/2013/Hume_2013_A*"
for sp in species.Species
    ismissing(sp) && continue
    run(`evince --find=$(sp) "/home/raf/zotero_library/Journal Article/undefined/2013/Hume_2013_A synopsis of the pre-human avifauna of the Mascarene Islands.pdf"`)
end


