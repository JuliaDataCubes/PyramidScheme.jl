module PyramidSchemeTylerExt
import Tyler
import PyramidScheme: PyramidProvider, selectlevel_tile
import MapTiles: Tile, extent
import Colors: RGB
import ColorSchemes: colorschemes
import DimensionalData: DimArray, Near
import TileProviders

function Tyler.fetch_tile(p::PyramidProvider, tile::Tile)
    ext = extent(tile, wgs84)
    data = try
        selectlevel_tile(p.p, ext, target_imsize=(256, 256))
    catch e
        println("Error getting extent $ext")
        println(e)
        return fill(RGB{Float64}(1.0, 1.0, 1.0), 256, 256)
    end
    if isempty(data)
        return fill(RGB{Float64}(1.0, 1.0, 1.0), 256, 256)
    end
    ar = DimArray(data.data, data.axes)
    cs = colorschemes[p.colorscheme]
    inner_lon_step = (ext.X[2] - ext.X[1]) / 256
    inner_lat_step = (ext.Y[2] - ext.Y[1]) / 256
    inner_lons = range(ext.X..., length=257) .+ inner_lon_step / 2
    inner_lats = range(ext.Y..., length=257) .+ inner_lat_step / 2
    map(CartesianIndices((256:-1:1, 1:256))) do I
        ilat, ilon = I.I
        v = ar[lon=Near(inner_lons[ilon]), lat=Near(inner_lats[ilat])] |> only
        if isnan(v) || ismissing(v) || isinf(v)
            p.nodatacolor
        else
            cs[(v-p.data_min)/(p.data_max-p.data_min)]
        end
    end
end
TileProviders.max_zoom(p::PyramidProvider) = p.max_zoom
TileProviders.min_zoom(p::PyramidProvider) = p.min_zoom

end