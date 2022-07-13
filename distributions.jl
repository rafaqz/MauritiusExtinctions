using Shapefile, Tables
using Plots
tbl = Shapefile.Table("/home/raf/PhD/Mauritius/Data/Distributions/REPTILES/REPTILES.shp") |> DataFrame

Phelsuma cepediana
x = filter(x -> occursin("bojeri", x.binomial, ), tbl)
plot(dems.mus)
plot!(x[1, :geometry])
