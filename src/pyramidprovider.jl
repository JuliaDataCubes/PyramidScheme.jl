import TileProviders
import Colors: RGB
import HTTP
import FileIO: save, Stream, @format_str
using MapTiles: WGS84, Tile
using ColorSchemes: colorschemes

struct PyramidProvider{P<:Pyramid{T} where T} <: TileProviders.AbstractProvider
    p::P
    min_zoom::Int
    max_zoom::Int
    data_min::Float64
    data_max::Float64
    colorscheme::Symbol
    nodatacolor::RGB{Float64}
end
PyramidProvider(p::Pyramid, data_min, data_max; min_zoom=0, max_zoom=nlevels(p), colorscheme=:viridis, nodatacolor=RGB{Float64}(1.0, 1.0, 1.0)) =
    PyramidProvider(p, min_zoom, max_zoom, data_min, data_max, colorscheme, nodatacolor)

function selectlevel_tile(pyramid, ext; target_imsize=(256, 256))
    pyrext = extent(pyramid)
    basepixels = map(pyrext, ext, size(pyramid)) do bbpyr, bbext, spyr
        pyrspan = bbpyr[2] - bbpyr[1]
        imsize = bbext[2] - bbext[1]
        imsize / pyrspan * spyr
    end
    dimlevels = log2.(basepixels ./ target_imsize)
    minlevel = minimum(dimlevels)
    n_agg = min(max(floor(Int, minlevel), 0), nlevels(pyramid))
    @debug "Selected level $n_agg"
    extcorrectnames = Extent{keys(pyrext)}(tuple(bounds(ext)...))
    levels(pyramid)[n_agg][extcorrectnames]
end

function provider_request_handler(p::PyramidProvider)
    request -> begin
        k = request.target
        #@show request.target
        k = lstrip(k, '/')
        m = match(r"(\d+)/(\d+)/(\d+).png", k)
        if m === nothing
            return HTTP.Response(404, "Error: Malformed url")
        end
        z, x, y = tryparse.(Int, m.captures)
        any(isnothing, (x, y, z)) && return HTTP.Response(404, "Error: Malformed url")
        data = fetch_tile(p, Tile(x, y, z))
        buf = IOBuffer()
        save(Stream{format"PNG"}(buf), data)
        return HTTP.Response(200, take!(buf))
    end
end

function fetch_tile(p::PyramidProvider, tile::Tile)
    ext = extent(tile, WGS84())
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
    ar = DD.DimArray(data.data, data.axes)
    cs = colorschemes[p.colorscheme]
    inner_lon_step = (ext.X[2] - ext.X[1]) / 256
    inner_lat_step = (ext.Y[2] - ext.Y[1]) / 256
    inner_lons = range(ext.X..., length=257) .+ inner_lon_step / 2
    inner_lats = range(ext.Y..., length=257) .+ inner_lat_step / 2
    map(CartesianIndices((256:-1:1, 1:256))) do I
        ilat, ilon = I.I
        v = ar[lon=DD.Near(inner_lons[ilon]), lat=DD.Near(inner_lats[ilat])] |> only
        if isnan(v) || ismissing(v) || isinf(v)
            p.nodatacolor
        else
            cs[(v-p.data_min)/(p.data_max-p.data_min)]
        end
    end
end
HTTP.serve(p::PyramidProvider, args...; kwargs...) = HTTP.serve(provider_request_handler(p), args...; kwargs...)
HTTP.serve!(p::PyramidProvider, args...; kwargs...) = HTTP.serve!(provider_request_handler(p), args...; kwargs...)
