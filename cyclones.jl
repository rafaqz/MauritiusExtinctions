using Raster
using CSV
using Shapefile
using TableView, Blink


cyclones = CSV.File("/home/raf/PhD/Mauritius/ibtracs.since1980.list.v04r00.csv"; skipto=3) |> DataFrame
cyclones
w = Blink.Window()
body!(w, showtable(cyclones))
names(cyclones)
wind_vars = names(cyclones)[occursin.(Ref("WIND"), names(cyclones))]
collect(skipmissing(cyclones.WMO_WIND))

world = Shapefile.Handle("/home/raf/PhD/Mauritius/world-administrative-boundaries/world-administrative-boundaries.shp")



function fill_cyclones!(A, cyclones)
    hurricane = 0
    prevlon = -1000.0
    prevlat = -1000.0
    new_hurricane = true
    fillval = Ref(0.0)
    fill(x) = max(x, fillval[])
    for row in eachrow(cyclones)
        windvals = [
            row[:WMO_WIND],
            row[:USA_WIND],
            row[:TOKYO_WIND],
            row[:CMA_WIND],
            row[:HKO_WIND],
            row[:NEWDELHI_WIND],
            row[:REUNION_WIND],
            row[:BOM_WIND],
            row[:NADI_WIND],
            row[:WELLINGTON_WIND],
            row[:DS824_WIND],
            row[:TD9636_WIND],
            row[:TD9635_WIND],
            row[:NEUMANN_WIND],
            row[:MLC_WIND],
        ]
        i = findfirst(x -> x != " ", windvals)
        isnothing(i) && continue
        lon = row.LON
        lat = row.LAT
        if hurricane != row.SID 
            hurricane = row.SID
        else
            line = (start=(x=prevlon, y=prevlat), stop=(x=lon, y=lat))
            length = sqrt((line.start.x - line.stop.x)^2 + (line.start.y - line.stop.y)^2)
            # if length > 3
                # @show hurricane line
            # else
                fillval[] = Base.parse(Float64, windvals[i])
                Rasters._fill_line!(A, line, fill, (X(), Y()))
            # end
        end
        # cycloneraster[X=Contains(lon), Y=Contains(lat)] = Base.parse(Float64, windval)
        # @show lon row.LON row.LAT 
        prevlon = lon
        prevlat = lat
    end
end

cyclone_long = zeros(X(Projected(-360:1:359.5; sampling=Intervals(Start()), crs=EPSG(4326))),
                     Y(Projected(-90:1:89.5; sampling=Intervals(Start()), crs=EPSG(4326))))

fill_cyclones!(cyclone_long, cyclones)

cycloneraster = cyclone_long[X(180:540)]
view(cycloneraster, X(1:180)) .= max.(view(cyclone_long, X(541:720)), view(cycloneraster, X(1:180)))
view(cycloneraster, X(181:360)) .= max.(view(cyclone_long, X(1:180)), view(cycloneraster, X(181:360)))

Plots.plot(cycloneraster .^ 2)
Plots.plot!(world; fill=nothing)

rasterize!(cycloneraster, world.shapes; fill=1, shape=:line)


