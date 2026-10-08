"""
PyramidScheme is a package to generate and work with Pyramids. 
A pyramid is a data structure that holds a two dimensional base array and corresponding aggregations,
which can be used for interactive plotting and analysis.
The main entry point of the package is the `Pyramid` type for the use of already available pyramids
and the `buildpyramids` function to generate the aggregation layers for an existing dataset.
"""
module PyramidScheme

using DiskArrayEngine: DiskArrayEngine, GMDWop, InputArray, LocalRunner, RegularWindows
using DiskArrayEngine: create_outwindows, engine
using DiskArrays: DiskArrays
using YAXArrayBase: YAXArrayBase as YAB
using YAXArrays: Cube, YAXArray, to_dataset, savedataset, savecube, setchunks, open_dataset
using ZarrCore: ZarrCore, zcreate, zopen, zgroup
using ZarrZstd: ZarrZstd
# Loaded so that Blosc compressed stores, e.g. older pyramids, can be read
using ZarrBlosc: ZarrBlosc
using DimensionalData: DimensionalData as DD
using DimensionalData.Dimensions: XDim, YDim
using Extents: Extent, extent, intersects
using FillArrays: Zeros
using Proj: Proj
using OffsetArrays: OffsetArray

using Statistics: mean, median

export Pyramid, buildpyramids, ConventionStyle, CellCenterStyle

"""
    Pyramid
A Pyramid will act as a DimArray with the data of the highest resolution,
but subsetting will return another Pyramid.
"""
struct Pyramid{T,N,D,A,B<:DD.AbstractDimArray{T,N,D,A},L, Me} <: DD.AbstractDimArray{T,N,D,A}
    "Highest resolution data of the pyramid"
    base::B
    "Aggregation layers of the `base` layer."
    levels::L
    "Metadata that describes the aggregation step if available"
    metadata::Me
end

"""
    Pyramid(data::AbstractDimArray; resampling_method=mean ∘ skipmissing, nlevels=nothing, factors=nothing, edges=:keep, recursive=true, tilesize=256, spatial_dims=SpatialDim)

Compute the pyramid levels of `data` in memory by aggregating windows of the spatial dimensions with `resampling_method`.
The levels are given by `nlevels` or `factors`, see [`pyramid_factors`](@ref).
Without either, enough levels are computed for the coarsest one to fit into `tilesize`.
`edges` decides whether partial windows at the end of a dimension are kept as their own cell (`:keep`)
or left out (`:drop`).
"""
function Pyramid(data::DD.AbstractDimArray; resampling_method= mean ∘ skipmissing, kwargs...)
    pyrdata, pyraxs = getpyramids(resampling_method, data; kwargs...)
    levels = DD.DimArray.(pyrdata, pyraxs)
    meta = Dict(deepcopy(DD.metadata(data)))
    push!(meta, "resampling_method" => resampling_method_name(resampling_method))
    Pyramid(data, levels, meta)
end

Pyramid(path::AbstractString) = Pyramid(path, YAB.backendfrompath(path)(path))
function Pyramid(path::AbstractString, backend)
    #This should rather be solved via dispatch, but this is not working because of Requires in YAXArrayBase.
    if backend isa YAB.backendlist[:zarr]
        _pyramid_zarr(path)
    elseif backend isa YAB.backendlist[:gdal]
        _pyramid_gdal(path)
    else
        throw(ArgumentError("""
        Loading is only supported for Zarr and GDAL Datasets got $backend.
        If you want to use GDAL you first have to load ArchGDAL.jl    
        """))
    end
end
function _pyramid_gdal end

function _pyramid_zarr(path)
    g = zopen(path)
    allkeys = collect(keys(g.groups))
    base = Cube(path)
    levavail = extrema(parse.(Int,allkeys[contains.(allkeys, r"\d")]))
    clevels = [Cube(open_dataset(g[string(l)])) for l in 1:last(levavail)]
    Pyramid(base, clevels, DD.metadata(base))
end
# refdims
# name

"""
    levels(pyramid::Pyramid)
Return all levels of the `pyramid`. These are order from the base to the coarsest aggregation.
This is an OffsetArray starting at zero for the base of the pyramid.
"""
levels(pyramid::Pyramid) = OffsetArray([pyramid.base, pyramid.levels...], 0:length(pyramid.levels))

levels(pyramid::Pyramid, i::Integer) = i==0 ? pyramid.base : pyramid.levels[i]
"""
    nlevels(pyramid)
Return the number of levels of the `pyramid`
"""
nlevels(pyramid::Pyramid) = length(levels(pyramid)) - 1
Base.parent(pyramid::Pyramid) = pyramid.base
Base.size(pyramid::Pyramid) = size(parent(pyramid))
function Base.:(==)(pyrA::Pyramid, pyrB::Pyramid)
    nlevela = nlevels(pyrA)
    nlevelb = nlevels(pyrB)
    nlevela != nlevelb && return false
    for i in nlevela:-1:1
        isequal(pyrA.levels[i], pyrB.levels[i]) || return false
    end
    return isequal(pyrA.base, pyrB.base)
end

Base.isequal(pyrA::Pyramid, pyrB::Pyramid) = pyrA == pyrB

#   all(isequal.(reverse(PyramidScheme.levels(pyrA)), reverse(PyramidScheme.levels(pyrB))))

DD.name(pyramid::Pyramid) = DD.name(parent(pyramid))
DD.refdims(pyramid::Pyramid) = DD.refdims(parent(pyramid))
DD.dims(pyramid::Pyramid) = DD.dims(parent(pyramid))
DD.metadata(pyramid::Pyramid) = pyramid.metadata

function DD.modify(f, pyr::Pyramid)
    pbase = DD.modify(f, pyr.base)
    plevels = map(pyr.levels) do level
        DD.modify(f, level)
    end
    Pyramid(pbase, plevels, pyr.metadata)
end
Base.read(pyr::Pyramid) = DD.modify(Array, pyr)
@inline function DD.rebuild(A::Pyramid, data, dims::Tuple=DD.dims(A), refdims=DD.refdims(A), name=DD.name(A))
    Pyramid(DD.rebuild(parent(A), data, dims, refdims, name, nothing), A.levels, A.metadata)
end

@inline function DD.rebuild(A::Pyramid; kw...)
    newbase = DD.rebuild(parent(A);dims=DD.dims(A), kw...)
    Pyramid(newbase, A.levels, DD.metadata(A))
end



@inline function DD.rebuildsliced(f::Function, A::Pyramid, data::AbstractArray, I::Tuple, name=DD.name(A))
    newbase = DD.rebuild(parent(A), parent(data), DD.slicedims(f, A, I)..., name)
    newlevels = map(A.levels) do level
        factors = map(level_factor, DD.dims(parent(A)), DD.dims(level))
        Ilevel = map(levelindex, factors, I, size(level))
        leveldata = f(parent(level), Ilevel...)
        DD.rebuild(level, leveldata, DD.slicedims(f, level, Ilevel)..., name)
    end
    Pyramid(newbase, newlevels, DD.metadata(A))
end

function Base.map(f, A1::Pyramid, As::Pyramid...)
    newbase = map(f, parent(A1), parent.(As)...)
    newlevels = [map(f, levels(A1,i), levels.(As, i)...) for i in 1:nlevels(A1)]
    Pyramid(newbase, newlevels, DD.metadata(A1)) # This should handle metadata better.
end

function DD.show_after(io::IO, mime, A::Pyramid)
    blockwidth = get(io, :blockwidth, 0)
    DD.print_block_separator(io, "pyramidlevels", blockwidth)
    println(io)
    println(io, "  Number of levels: $(nlevels(A)) ")
    for l in levels(A)
        println(io, "   ", size(l))
    end
    DD.print_block_close(io, blockwidth)
end
"""
Internal function

# Extended help

    aggregate_by_factor(xout, x, f, fac)


Aggregate windows of `fac` by `fac` cells of the data `x` using the function `f` and store the results in `xout`.
Windows at the end of `x` may be smaller.
This function is used as the inner function for the DiskArrayEngine call.
"""
function aggregate_by_factor(xout,x,f,fac)
    for j in axes(xout,2)
        for i in axes(xout,1)
            xview = ((i-1)*fac+1):min(size(x,1),(i*fac))
            yview = ((j-1)*fac+1):min(size(x,2),(j*fac))
            xout[i,j] = f(view(x,xview,yview))
        end
    end
end


"""
    all_pyramids!(xout, x, recursive, f, factors)
Compute all tiles of the pyramids for data `x` and store it in `xout`.
Uses function `f` as a aggregating function.
`factors` are the cumulative downsampling factors of the levels, starting with 1 for `x`.
`recursive` indicates whether higher tiles are computed from lower tiles or directly from the original data.
This is an optimization which for functions like median might lead to misleading results.
"""
function all_pyramids!(xout,x,recursive,f,factors)
    xnow = x
    for (i, xcur) in enumerate(xout)
        fac = recursive ? factors[i+1] ÷ factors[i] : factors[i+1]
        aggregate_by_factor(xcur,xnow,f,fac)
        if recursive
            xnow = xcur
        end
    end
end

"""
    gen_pyr(x...; recursive=true)
Generate the pyramid levels. `x` holds the output arrays of all levels, the input array,
the aggregation function and the cumulative downsampling factors, in this order.
"""
function gen_pyr(x...;recursive=true)
    arrays = x[1:end-2]
    f, factors = x[end-1], x[end]
    all_pyramids!(Base.front(arrays),last(arrays),recursive,f,factors)
end


"""
    fill_pyramids(data, outputs, func, recursive, spatial_dims, factors)
Fill the pyramids generated from the `data` with the aggregation function `func` into the list `outputs`.
`factors` are the cumulative downsampling factors, starting with 1 for `data`,
and `outputs` must have the sizes of the levels with partial windows kept, see [`level_length`](@ref).
`recursive` indicates whether higher tiles are computed from lower tiles or directly from the original data.
This is an optimization which for functions like median might lead to misleading results.
"""
function fill_pyramids(data, outputs,func,recursive,spatial_dims,factors;runner=LocalRunner, kwargs...)
    n_level = length(outputs)
    pixel_base_size = last(factors)
    pyramid_sizes = size.(outputs)
    tmp_sizes = [pixel_base_size ÷ f for f in factors[2:end]]

    inputwindows = Base.OneTo.(size(data))
    outputwindows = Base.OneTo.(size(data))
    for d in spatial_dims
        inputwindows = Base.setindex(inputwindows,RegularWindows(1,size(data,d),window=pixel_base_size),d)
    end

    ia  = InputArray(data, windows = inputwindows)


    oa = ntuple(n_level) do i
        for d in spatial_dims
            outputwindows = Base.setindex(outputwindows,RegularWindows(1,pyramid_sizes[i][d],window=tmp_sizes[i]),d)
        end
        create_outwindows(pyramid_sizes[i],windows = outputwindows)
    end

    func = DiskArrayEngine.create_userfunction(gen_pyr,ntuple(_->eltype(first(outputs)),length(outputs));is_mutating=true,kwargs = (;recursive),args = (func, factors))


    op = GMDWop((ia,), oa, func)

    lr = DiskArrayEngine.optimize_loopranges(op,5e8,tol_low=0.2,tol_high=0.05,max_order=2)
    r = runner(op,lr,outputs;kwargs...)
    run(r)
end

"""
    ESALCMode(counts)

Struct for more efficient handling of aggregating categorical land cover data.
"""
struct ESALCMode
    counts::Vector{Vector{Int}}
end

ESALCMode() = ESALCMode([zeros(Int,256) for _ in 1:Threads.nthreads()])

function (f::ESALCMode)(x)
    cv = f.counts[Threads.threadid()]
    fill!(cv,0)
    for ix in x
        cv[Int(ix)+1] += 1
    end
    _,ind = findmax(cv)
    UInt8(ind-1)
end


"""
    arraywindows(s,w)

Construct a list of `RegularWindows` for the size list in `s` for windows `w`.
??
"""
function arraywindows(s,w)
    map(s) do l
        RegularWindows(1,l,window=w)
    end
end


"""
    compute_nlevels(data, tilesize=1024)

Compute the number of levels for the aggregation based on the size of `data`.
"""
compute_nlevels(data, tilesize=256) = max(0,ceil(Int,log2(maximum(size(data))/tilesize)))

"""
    gen_output(t,s)

Create output array of type `t` and size `s`
If the array is smaller than 10e6 it is created on disk and otherwise as a temporary zarr file.
"""
function gen_output(t,s; path=tempname())
    # This should be dispatching on the output type whether it is internal or external
    outsize = sizeof(t)*prod(s)
    if outsize > 100e6
        # This should be zgroup instead of zcreate, could use savedataset(skelton=true)
        # Dummy dataset with FillArrays with the shape of the pyramidlevel
        zcreate(t,s...;path,chunks = (1024,1024),fill_value=zero(t))
    else
        zeros(t,s...)
    end
end
function DiskArrayEngine.engine(dimarr::DD.AbstractDimArray) 
    dataengine = engine(dimarr.data)
    DD.rebuild(dimarr, data=dataengine)
end
DiskArrayEngine.engine(pyr::Pyramid) = Pyramid(engine(parent(pyr)), engine.(pyr.levels), DD.metadata(pyr))
"""
    output_arrays(pyramid_sizes)

Create the output arrays for the given `pyramid_sizes`
"""
output_arrays(pyramid_sizes, T) = [gen_output(T,p) for p in pyramid_sizes]

"""
    SpatialDim
Union of Dimensions which are assumed to be in space and are therefore used in the pyramid building. 
"""
SpatialDim = (DD.Dimensions.XDim, DD.Dimensions.YDim)

"""
    getpyramids(reducefunc, ras; recursive=true, tilesize=256, nlevels=nothing, factors=nothing, edges=:keep, spatial_dims=SpatialDim)
Compute the data of the pyramids of a given data cube `ras`.
This returns the data of the pyramids and the dimensions of the levels.

The levels are given by `nlevels` or `factors`, see [`pyramid_factors`](@ref).
Without either, enough levels are computed for the coarsest one to fit into `tilesize`.
`edges` decides whether partial windows at the end of a dimension are kept as their own cell
(`:keep`) or left out (`:drop`).
"""
function getpyramids(reducefunc, ras;recursive=true, tilesize=256, nlevels=nothing, factors=nothing, edges=:keep, spatial_dims=SpatialDim)
    factors = pyramid_factors(; nlevels, factors, default_nlevels=compute_nlevels(ras, tilesize))
    if length(factors) == 1
        @info "Array is smaller than the tilesize no pyramidlevels are computed"
        return [], []
    end
    pyramid_axes = [level_dims(DD.dims(ras), f; edges, spatial_dims) for f in factors[2:end]]
    # Levels are computed with partial windows and trimmed afterwards,
    # the cells of partial windows do not contribute to coarser levels inside the trimmed sizes.
    keep_sizes = [length.(level_dims(DD.dims(ras), f; edges=:keep, spatial_dims)) for f in factors[2:end]]
    outtype = Base.infer_return_type(reducefunc, (Matrix{eltype(ras)},))
    #outtype = Missing <: eltype(ras) ? Union{Missing, outtypefunc} : outtypefunc
    outmin = output_arrays(keep_sizes, outtype)
    ispatial_dims = DD.dimnum(DD.dims(ras),spatial_dims)

    fill_pyramids(ras,outmin,reducefunc,recursive, ispatial_dims, factors; threaded=true)

    leveldata = map(outmin, pyramid_axes) do out, axes
        size(out) == length.(axes) ? out : view(out, Base.OneTo.(length.(axes))...)
    end
    leveldata, pyramid_axes
end

"""
    selectlevel(pyramids, ext, resolution; target_imsize=(1024,512)
Internal function to select the raster data that should be plotted on screen. 
`pyramids` is a Vector of Raster objects with increasing coarsity. 
`ext` is the extent of the zoomed in area
`target_imsize` is the target size of the output data that should be plotted.
"""
function selectlevel(pyramid, ext;target_imsize=(1024, 512))
    pyrext = extent(pyramid)
    intersects(pyrext, ext) || return zero(pyramid.levels[end])
    basepixels = map(keys(pyrext)) do bb
        pyrspan = pyrext[bb][2] - pyrext[bb][1]
        imsize = ext[bb][2] - ext[bb][1]
        imsize / pyrspan * size(pyramid, bb)
    end
    dimlevels = log2.(basepixels ./ target_imsize)
    minlevel = maximum(dimlevels)
    n_agg = min(max(ceil(Int,minlevel),0),nlevels(pyramid))
    @debug "Selected level $n_agg"
    levels(pyramid)[n_agg][ext]
end




xkey(keyext) = DD.name(DD.dims(keyext, XDim))
ykey(keyext) = DD.name(DD.dims(keyext, YDim))
#TODO write test, move to utils.jl

"""
    switchkeys(dataext, keyext)
Internal function

### Extended help
    Return an Extent with the limits from `dataext` and the keys of `keyext`.
    We assume that dataext has keys X, and Y and the keys of keyext are XDim and YDim from DimensionalData
"""
function switchkeys(dataext, keyext)
    xk = xkey(keyext)
    yk = ykey(keyext)
    nt = NamedTuple{(xk, yk)}((dataext.X, dataext.Y))
    Extent(nt)
end

"""
    trans_bounds(trans::Proj.Transformation, bbox::Extent, densify_pts::Integer)
Compute the projection of `bbox` with the transformation `trans`. 
This is used to project the data on the fly to another transformation.
"""
function trans_bounds(
    trans::Proj.Transformation,
    bbox::Extent,
    densify_pts::Integer = 21,
)::Extent
    xlims, ylims = Proj.bounds(trans, bbox.X, bbox.Y; densify_pts)
    return Extent(X = xlims, Y = ylims)
end

function Base.cat(A1::Pyramid, As::Pyramid...;dims)
    catlevels = [cat(A1.levels[i], levels.(As, i)...; dims) for i in eachindex(A1.levels)]
    catbase = cat(parent(A1), parent.(As)...; dims)
    Pyramid(catbase, catlevels, merge(DD.metadata(A1), DD.metadata.(As)...))
end



include("conventions.jl")
include("levels.jl")
include("writer.jl")
include("broadcast.jl")
end
