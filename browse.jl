using DataFrames, CSV

# Google search all species, opening the next when the browser is closed.
for row in eachrow(s)
    # ismissing(row["Species"]) && continue
    # ismissing(row["Common_name"]) && continue
    row.Origin == "Endemic" || continue

    # placenamed = false
    # for place in ("RéuReunion", "Mauritius", "Rodrigues", "Mascarene")
    #     placenamed = occursin(place, row["Common name"]) 
    #     @show row["Common name"] place placenamed
    #     placenamed && break
    # end
    # placenamed && continue
    search = replace(row["Species"], " " => "+")
    # search = row["Species"]
    # run(`chromium https\://www.google.com/search\?q=\"$(search)\"`)
    run(`chromium https://www.iucnredlist.org/search\?query=$search\&searchType=species`)
    # run(`chromium https\://www.google.com/search\?tbm=isch\&q=$(search)`)
end

for row in Tables.rows(mascarine_species)
    sp = row.Species
    row.Origin === "Alien" && continue
    (!ismissing(row.Reunion) && row.Reunion == "TRUE") || continue
    ismissing(sp) && continue
    search = replace(sp, " " => "+")
    clipboard(row.Common_name)
    run(`chromium https\://www.google.com/search\?q=$("search")\&tbm=isch\&tbs=isz:cl`)
    # run(`chromium https\://www.google.com/search\?tbm=isch\&q=$(search)`)
end

# Check in Hume
humepath = "/home/raf/zotero_library/Journal Article/undefined/2013/Hume_2013_A*"
for sp in species.Species
    ismissing(sp) && continue
    run(`evince --find=$(sp) "/home/raf/zotero_library/Journal Article/undefined/2013/Hume_2013_A synopsis of the pre-human avifauna of the Mascarene Islands.pdf"`)
end

plantpath = "/home/raf/Data/Extinction/Plants/Mascarene_extinct_plants.csv"
extinct_plants = CSV.read(plantpath, DataFrame)

for row in Tables.rows(extinct_plants)
    sp = row["Accepted.binomial"]
    search = replace(sp, "_" => "+")
    # run(`chromium https\://www.google.com/search\?q=$(search)\&tbm=isch\&tbs=isz:cl`)
    run(`chromium https\://www.google.com/search\?q=$(search)`)
end
