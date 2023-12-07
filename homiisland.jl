homiisland = Raster("/home/raf/PhD/Mascarenes/Data/Generated/Landcover/mus_landcover.tif")
homiisland_resampled = resample(homiisland; to=dems.mus)
Rasters.rplot(homiisland_resampled .== 8)
write("/home/raf/PhD/Mascarenes/Data/Generated/Landcover/mus_landcover_resampled.tif")
