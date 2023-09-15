homiisland = Raster("/home/raf/PhD/Mascarenes/Data/Generated/Landcover/mus_landcover.tif")
homiisland2 = resample(homiisland; to=dems.mus)
Rasters.rplot(homiisland2 .== 8)
write("/home/raf/PhD/Mascarenes/Data/Generated/Landcover/mus_landcover_2.tif")
