using CSV
using DataFrames
using TableView
using Blink

pan_theria = CSV.read("/home/raf/PhD/Traits/PanTHERIA/ECOL_90_184/PanTHERIA_1-0_WR05_Aug2008.txt", DataFrame; missingstring=["-999.00", "-999"])

p_niger = subset(pan_theria, :MSW05_Genus => ByRow(==("Pteropus")), :MSW05_Species => ByRow(==("niger")))[1, :]
filter(x -> !ismissing(x[2]), names(p_niger) .=> values(p_niger))

elton_birds = CSV.read("/home/raf/PhD/Traits/EltonTraits/BirdFuncDat.txt", DataFrame)
elton_mammals = CSV.read("/home/raf/PhD/Traits/EltonTraits/MamFuncDat.txt", DataFrame)

lizards = CSV.File("/home/raf/PhD/Traits/Lizards/Appendix S1 - Lizard data version 1.0.csv") |> DataFrame

w = Blink.Window()
interface = TableView.showtable(lizards[1:1000, :]; height=1000);
Blink.body!(w, interface)
