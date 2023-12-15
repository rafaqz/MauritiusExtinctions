using JSON3, MapRasterization, GeoInterface, Rasters, FileIO, ImageIO, DataFrames, CSV, GeoJSON
using DimensionalData.LookupArrays

function open_output(T, filename::String)
    json_path = splitext(filename)[1] * ".json"
    return isfile(json_path) ? JSON3.read(read(json_path), T) : nothing
end

open_warp_points(x::NamedTuple) = open_warp_points(x.filename)
function open_warp_points(filename::String)
    csv_path = splitext(filename)[1] * ".csv"
    return isfile(csv_path) ? CSV.read(csv_path, DataFrame) : nothing
end

function warp_to_raster(img_path::String, template::Raster;
    object_type=MapRasterization.MapSelection, edit=false, save=true, kw...
)
    img = load_image(img_path)
    csv_path = splitext(img_path)[1] * ".csv"
    points = isfile(csv_path) ? open_warp_points(img_path) : nothing
    if edit || !isfile(csv_path)
        warp_points = if isnothing(points)
            MapRasterization.click_warp(img;
                template=reverse(template; dims=Y()), kw...
            )
        else
            :x_known in names(points) && rename!(points, [:x_known => :x_a, :y_known => :y_a, :x_unknown => :x_b, :y_unknown => :y_b])
            # if :x_a in names(points)
                MapRasterization.click_warp(img;
                    template=reverse(template; dims=Y()), points, kw...
                )
            # else
                # MapRasterization.click_warp(Float64.(Gray.(img));
                    # template=reverse(template; dims=Y()), missingval=missingval(template),
                # )
            # end
        end
        if save
            df = DataFrame(warp_points)
            CSV.write(csv_path, df)
        end
    end
    if save && isfile(splitext(img_path)[1] * ".json")
        df = CSV.read(csv_path, DataFrame)
        output = open_output(object_type, img_path)
        poly = 1
        if object_type <: MapRasterization.MapSelection
            out = Int.(reshape(output.output, size(img)))
            warper = MapRasterization.Warper(df, template, poly)
            rs = MapRasterization.warp(warper, out; missingval=0)
            mask!(rs; with=template)
            raster_path = splitext(img_path)[1] * ".tif"
            write(raster_path, rs; force=true)
            return rs
        elseif object_type <: MapRasterization.CategoryShapes
            warper = MapRasterization.Warper(df, template, poly)
            warped_geoms = map(output.shapes) do sh
                geoms = Polygon.(sh)
                warped = MapRasterization.warp(warper, geoms)
                map(warped) do g
                    collect(GeoInterface.getpoint(g))
                end
            end
            warped = MapRasterization.CategoryShapes{Polygon}(warped_geoms, output.category_names)
            # GeoJSON.write(splitext(img_path)[1] * "_warped.geojson", warped)
            return warped
        end
    end
end

function choose_categories(img_path::String;
    save=true, restart=false,
    output=restart ? nothing : open_output(MapRasterization.MapSelection, img_path),
)
    img = load_image(img_path)
    if isnothing(output)
        cs = MapRasterization.CategorySelector(img)
    else
        cs = MapRasterization.CategorySelector(img, output)
    end
    # if save
    #     output = MapRasterization.MapSelection(cs)
    #     json_path = splitext(img_path)[1] * ".json"
    #     if isfile(json_path)
    #         backup_path = splitext(img_path)[1] * "_backup.json"
    #         cp(json_path, backup_path; force=true)
    #     end
    #     write(json_path, JSON3.write(output))
    # end
    return cs
end

load_image(img_path::String) = RGB{Float64}.(load(img_path) |> rotr90)
