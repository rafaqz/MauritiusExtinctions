using Rasters, GLMakie
includet("raster_common.jl")

# Build auxiliary rasters
lc_predictions_paths = (
    mus="$outputdir/lc_predictions_mus.nc",
    reu="$outputdir/lc_predictions_reu.nc",
    rod="$outputdir/lc_predictions_rod.nc",
)
mus_veg_path = "/home/raf/PhD/Mascarenes/Data/Selected/Mauritius/Undigitised/page33_mauritius_vegetation_colored.tif"
reu_veg_path = "/home/raf/PhD/Mascarenes/Data/Dominique/Vegetation_Rasters/pastveg3.tif"
using DBFTables
original_veg = (;
    mus=reorder(replace_missing(Raster(mus_veg_path), 0), lc_predictions.mus),
    reu=reorder(resample(replace_missing(Raster(reu_veg_path), 0); to=masks.reu), lc_predictions.reu),
)
# netcdf has the annoying center locus for time
lc_predictions = map(lc_predictions_paths) do path
    lc_predictions = RasterStack(path) |>
        x -> rebuild(Rasters.modify(BitArray, x); missingval=false) |>
        x -> Rasters.set(x, Ti => Int.(maybeshiftlocus(Start(), dims(x, Ti), )))
end
original_veg.mus
lc_predictions.mus

island_veg_change = map(lc_predictions[keys(original_veg)], original_veg) do p, v
    rebuild(UInt8.(broadcast_dims(*, p.native, v)); missingval=0)
end
p = Rasters.rplot(island_veg_change.mus[Ti=At(1700:2018)]; colormap=:viridis)
save("images/mus_original_veg.png", p)


reu_veg_df = DBFTables.Table("/home/raf/PhD/Mascarenes/Data/Dominique/Vegetation_Rasters/pastveg3.tif.vat.dbf") |> DataFrame
island_habitat_names = (;
    mus=["semi-dry_evergreen_forest", "open_dry_palm-rich_woodland", "wet_forest", "pandanus_swamp", "mossy_rainforest", "mangrove", "wetland vegetation"],
    reu=reu_veg_df.HABITAT,
)
nhabitats = (mus=7, reu=20)
habitat_sums = map(island_veg_change, island_habitat_names, nhabitats) do veg_change, habitat_names, nhabitat
    As = map(1:nhabitat) do habitat
        dropdims(sum(==(habitat), veg_change; dims=(X, Y)); dims=(X, Y))
    end
    cat(As...; dims=Dim{:habitat}(habitat_names))
end
k = :reu
k = :mus
cum = cumsum(habitat_sums[k]; dims=2)
x = lookup(habitat_sums[k], Ti)
fig = Figure()
ax = Axis(fig[1, 1])
for i in nhabitats[k]:-1:1
    y = parent(cum[habitat=i])
    Makie.lines!(x, y; color=:black)
    band!(x, fill(0, length(x)), y; label = "Label")
end
fig[1, 2] = Legend(fig, ax, habitat_names)


