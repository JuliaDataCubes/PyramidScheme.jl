import DiskArrays
import TileProviders: AbstractProvider
import Colors: RGB, RGBA
import HTTP
import FileIO
import MapTiles: Tile, extent, web_mercator
import DimensionalData: Dim, X, Y

"DiskArray representing maptiles as 2d or 3d arrays with bands"
struct MapTileDiskArray{T,N,P<:AbstractProvider} <: DiskArrays.ChunkTiledDiskArray{T,N}
    prov::P
    zoom::Int
    tilesize::Tuple{Int,Int}
    nband::Int
end
function MapTileDiskArray(prov, zoom, mode=:band)
    testtile = load_data(prov, zoom, 0, 0)
    et = eltype(testtile)
    nband = if et <: RGB
        3
    elseif et <: RGBA
        4
    else
        error("Unknown color type $et")
    end
    if mode === :band
        return MapTileDiskArray{rgbeltype(et),3,typeof(prov)}(prov, zoom, size(testtile, 1), nband)
    elseif mode === :rgb
        return MapTileDiskArray{et,2,typeof(prov)}(prov, zoom, size(testtile, 1), nband)
    else
        error("Unknown mode $mode")
    end
end
DiskArrays.eachchunk(a::MapTileDiskArray{<:Any,3}) = DiskArrays.GridChunks((a.nband, a.tilesize[1] * 2^a.zoom, a.tilesize[2] * 2^a.zoom), (a.nband, a.tilesize...))
DiskArrays.eachchunk(a::MapTileDiskArray{<:Any,2}) = DiskArrays.GridChunks((a.tilesize[1] * 2^a.zoom, a.tilesize[2] * 2^a.zoom), a.tilesize)

DiskArrays.haschunks(::MapTileDiskArray) = DiskArrays.Chunked()

rgbeltype(::Type{RGB{T}}) where {T} = T
rgbeltype(::Type{RGBA{T}}) where {T} = T
function Base.getindex(a::MapTileDiskArray{<:Any,3}, i::DiskArrays.ChunkIndex{3,DiskArrays.OneBasedChunks})
    _, y, x = i.I.I
    data = load_data(a.prov, a.zoom, x - 1, y - 1)
    T = rgbeltype(eltype(data))
    data = reinterpret(reshape, T, data)
    return DiskArrays.wrapchunk(data, DiskArrays.eachchunk(a)[i.I])
end

function Base.getindex(a::MapTileDiskArray{<:Any,2}, i::DiskArrays.ChunkIndex{2,DiskArrays.OneBasedChunks})
    y, x = i.I.I
    data = load_data(a.prov, a.zoom, x - 1, y - 1)
    return DiskArrays.wrapchunk(data, DiskArrays.eachchunk(a)[i.I])
end

function load_data(prov, zoom, x, y)
    url = TileProviders.geturl(prov, x, y, zoom)
    result = HTTP.get(url; retry=false, readtimeout=4, connect_timeout=4)
    if result.status > 300
        if result.status == 404
            return nothing
        else
            throw(ErrorException("HTTP error $(result.status)"))
        end
    else
        io = IOBuffer(result.body)
        format = FileIO.query(io)
        FileIO.load(format)
    end
end

function dimsfromzoomlevel(zoom, tilesize)
    t1 = Tile(1, 1, zoom)
    ntiles = 2^zoom
    npix = ntiles .* tilesize
    t2 = Tile(ntiles, ntiles, zoom)
    ex1 = extent(t1, web_mercator)
    ex2 = extent(t2, web_mercator)
    x1, x2 = first(ex1.X), last(ex2.X)
    y1, y2 = first(ex1.Y), last(ex2.Y)
    stepx = (x2 - x1) / npix[1]
    stepy = (y2 - y1) / npix[2]
    x = X(DD.Sampled(range(x1, x2 - stepx, length=npix[1]), sampling=DD.Intervals(DD.Start())))
    y = Y(DD.Sampled(range(y1, y2 - stepy, length=npix[2]), sampling=DD.Intervals(DD.Start())))
    return x, y
end
function provtoyax(prov, zoom, mode=:band)
    a = MapTileDiskArray(prov, zoom, mode)
    xdim, ydim = dimsfromzoomlevel(zoom, a.tilesize)
    if mode === :band
        coldim = if a.nband == 3
            Dim{:Band}(["Red", "Green", "Blue"])
        elseif a.nband == 4
            Dim{:Band}(["Red", "Green", "Blue", "Alpha"])
        end
        YAXArray((coldim, ydim, xdim), a)
    else
        YAXArray((ydim, xdim), a)
    end
end

function Pyramid(prov::AbstractProvider, mode=:band)
    maxzoom = get(prov.options, :max_zoom, 18)
    base = provtoyax(prov, maxzoom, mode)
    levels = [provtoyax(prov, zoom, mode) for zoom in (maxzoom-1):-1:0]
    return Pyramid(base, levels, prov.options)
end
