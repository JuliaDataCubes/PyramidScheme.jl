"""
PyramidScheme is a package to generate and work with Pyramids. 
A pyramid is a data structure that holds a two dimensional base array and corresponding aggregations,
which can be used for interactive plotting and analysis.
The main entry point of the package is the `Pyramid` type for the use of already available pyramids
and the `buildpyramids` function to generate the aggregation layers for an existing dataset.
"""
module PyramidScheme

using DiskArrayEngine: DiskArrayEngine, GMDWop, InputArray, LocalRunner, MovingWindow, RegularWindows
using DiskArrayEngine: create_outwindows, engine
using DiskArrays: DiskArrays
#using YAXArrays: savecube
using Zarr
using YAXArrayBase
const YAB = YAXArrayBase
using YAXArrays: Cube, YAXArray, to_dataset, savedataset, setchunks, open_dataset
using Zarr: zcreate, writeattrs
using DimensionalData: DimensionalData as DD
using DimensionalData.Dimensions: XDim, YDim
using Extents
using FillArrays: Zeros
using Proj
using Makie: Axis, Colorbar, DataAspect, Figure, FigureAxisPlot, Observable, on, heatmap!
import MakieCore: plot, plot!
using OffsetArrays

using Statistics



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

function Pyramid(data::DD.AbstractDimArray)
    pyrdata, pyraxs = getpyramids(mean ∘ skipmissing, data, recursive=false)
    levels = DD.DimArray.(pyrdata, pyraxs)
    meta = Dict(deepcopy(DD.metadata(data)))
    push!(meta, "resampling_method" => "mean_skipmissing")
    Pyramid(data, levels, meta)
end

Pyramid(path::AbstractString) = Pyramid(path, YAB.backendfrompath(path)(path))
function Pyramid(path::AbstractString, backend)
    #This should rather be solved via dispatch, but this is not working because of Requires in YAXArrayBase.
    if backend isa YAB.ZarrDataset
        _pyramid_zarr(path)
    elseif backend isa YAB.GDALDataset
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
    base = Cube(path)[Ti=1] # This getindex should be unnecessary and I should rather fix my data on disk
    levavail = extrema(parse.(Int,allkeys[contains.(allkeys, r"\d")]))
    clevels = [Cube(open_dataset(g[string(l)])) for l in 1:last(levavail)]
    Pyramid(base, clevels, Dict())
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

DD.name(pyramid::Pyramid) = DD.name(parent(pyramid))
DD.refdims(pyramid::Pyramid) = DD.refdims(parent(pyramid))
DD.dims(pyramid::Pyramid) = DD.dims(parent(pyramid))
DD.metadata(pyramid::Pyramid) = pyramid.metadata
@inline function DD.rebuild(A::Pyramid, data, dims::Tuple=dims(A), refdims=refdims(A), name=name(A))
    Pyramid(DD.rebuild(parent(A), data, dims, refdims, name, nothing), A.levels, A.metadata)
end

@inline DD.rebuild(A::Pyramid; kw...) = Pyramid(DD.rebuild(parent(A); kw...), A.levels, DD.metadata(A))


@inline function DD.rebuildsliced(f::Function, A::Pyramid, data::AbstractArray, I::Tuple, name=DD.name(A))
    newbase = DD.rebuild(parent(A), parent(data), DD.slicedims(f, A, I)..., name)
    newlevels = map(enumerate(A.levels)) do (z, level)
        Ilevel =  levelindex(z, I)
        leveldata = f(parent(level), Ilevel...)
        DD.rebuild(level, leveldata, DD.slicedims(f, level, Ilevel)..., name)
    end
    Pyramid(newbase, newlevels, DD.metadata(A))
end

function DD.show_after(io::IO, mime, A::Pyramid)
    blockwidth = get(io, :blockwidth, 0)
    DD.print_block_separator(io, "pyramidlevels", blockwidth)
    println(io, "  Number of levels: $(nlevels(A)) ")
    for l in levels(A)
        println(io, "   ", size(l))
    end
    DD.print_block_close(io, blockwidth)
end
"""
    levelindex(z, i)
Internal function.
# Extended help

Compute the index into the level `z` from the index i by shifting by a power of two on the Integer level
"""
levelindex(z, i::Integer) = (i - 1) >> z +1
# TODO: This should be also done generically for OrdinalRanges
levelindex(z, i::AbstractUnitRange) = levelindex(z, first(i)):levelindex(z, last(i))
levelindex(z, I::Tuple) = map(i -> levelindex(z, i), I)
"""
Internal function

# Extended help

    aggregate_by_factor(xout, x, f)


Aggregate the data `x` using the function `f` and store the results in `xout`.
This function is used as the inner function for the DiskArrayEngine call.
"""
function aggregate_by_factor(xout,x,f)
    fac = ceil(Int,size(x,1)/size(xout,1))
    for j in axes(xout,2)
        for i in axes(xout,1)
            xview = ((i-1)*fac+1):min(size(x,1),(i*fac))
            yview = ((j-1)*fac+1):min(size(x,2),(j*fac))
            xout[i,j] = f(view(x,xview,yview))
        end
    end
end


"""
    all_pyramids!(xout, x, recursive, f)
Compute all tiles of the pyramids for data `x` and store it in `xout`.
Uses function `f` as a aggregating function.
`recursive` indicates whether higher tiles are computed from lower tiles or directly from the original data. 
This is an optimization which for functions like median might lead to misleading results.
""" 
function all_pyramids!(xout,x,recursive,f)
    xnow = x
    for xcur in xout
        aggregate_by_factor(xcur,xnow,f)
        if recursive
            xnow = xcur
        end
    end
end

"""
    gen_pyr(x; recursive=true)
Generate the pyramid for the data `x`.
I do not understand what is `x`
"""
function gen_pyr(x...;recursive=true) 
    allx = Base.front(x)
    f = last(x)
    all_pyramids!(Base.front(allx),last(allx),recursive,f)
end


"""
    fill_pyramids(data, outputs, func, recursive)
Fill the pyramids generated from the `data` with the aggregation function `func` into the list `outputs`.
`recursive` indicates whether higher tiles are computed from lower tiles or directly from the original data. 
This is an optimization which for functions like median might lead to misleading results.
"""
function fill_pyramids(data, outputs,func,recursive;runner=LocalRunner, kwargs...)
    n_level = length(outputs)
    pixel_base_size = 2^n_level
    pyramid_sizes = size.(outputs)
    tmp_sizes = [ceil(Int,pixel_base_size / 2^i) for i in 1:n_level]

    ia  = InputArray(data, windows = arraywindows(size(data),pixel_base_size))

    oa = ntuple(i->create_outwindows(pyramid_sizes[i],windows = arraywindows(pyramid_sizes[i],tmp_sizes[i])),n_level)

    func = DiskArrayEngine.create_userfunction(gen_pyr,ntuple(_->eltype(first(outputs)),length(outputs));is_mutating=true,kwargs = (;recursive),args = (func,))

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

function agg_axis(d,n)
    DD.rebuild(d, LinRange(first(d), last(d), cld(length(d), n)))
end
"""
    gen_output(t,s)

Create output array of type `t` and size `s`
If the array is smaller than 10e6 it is created on disk and otherwise as a temporary zarr file.
"""
function gen_output(t,s; path=tempname())
    # This should be dispatching on the output type whether it is internal or external
    outsize = sizeof(t)*prod(s)
    z = Zeros(t, s...)
    if outsize > 100e6
        # This should be zgroup instead of zcreate, could use savedataset(skelton=true)
        # Dummy dataset with FillArrays with the shape of the pyramidlevel
        zcreate(t,s...,path=p,chunks = (1024,1024),fill_value=zero(t))
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
SpatialDim = Union{DD.Dimensions.XDim, DD.Dimensions.YDim}

"""
    buildpyramids(path; resampling_method=mean)
Build the pyramids for the zarr dataset at `path` and write the pyramid layers into the zarr folder.
The different scales are written according to the GeoZarr spec and a multiscales entry is added to the attributes of the zarr dataset.
The data is aggregated with the specified `resampling_method`. 
Keyword arguments are forwarded to the `fill_pyramids` function.
"""
function buildpyramids(path; resampling_method=mean, recursive=true, runner=LocalRunner, verbose=false)
    if YAB.backendfrompath(path) != YAB.ZarrDataset 
        throw(ArgumentError("$path  is not a Zarr dataset therefore we can't build the Pyramids inplace"))
    end

    # Should this be a dataset not a Cube?
    # Build a loop for all variables in a dataset?
    org = Cube(path)
    # We run the method once to derive the output type
    t = typeof(resampling_method(zeros(eltype(org), 2,2)))
    n_level = compute_nlevels(org)            
    input_axes = filter(x-> x isa SpatialDim,  DD.dims(org))
    if length(input_axes) != 2
        throw(ArgumentError("Expected two spatial dimensions got $input_axes"))
    end
    verbose && println("Constructing output arrays")
    outarrs = [output_zarr(n, input_axes, t, joinpath(path, string(n))) for n in 1:n_level]
    verbose && println("Start computation")
    fill_pyramids(org, outarrs, resampling_method, recursive;runner)
    pyraxs = [agg_axis.(input_axes, 2^n) for n in 1:n_level]
    pyrlevels = DD.DimArray.(outarrs, pyraxs)
    meta = Dict(deepcopy(DD.metadata(org)))
    push!(meta, "resampling_method" => string(resampling_method))
    multiscale = Dict{String, Any}()
    push!(multiscale, "datasets" => ["path"=> "", ["path" => string(i) for i in 1:n_level]...])
    push!(multiscale, "type" => "reduce")
    push!(meta, "multiscales" => [multiscale,])
    storage, basepath = Zarr.storefromstring(path)
    writeattrs(storage,basepath, meta)
    Pyramid(org, pyrlevels, meta)
    # Construct the TMS metadata from the info of the array
    # Save TMS to the .zattrs of the original data
end

# TODO: Decide on how we count levels rather from the bottom to the top or the other way around.
# I find bottom to top more intuitive but TMS specifies it differently

"""
    output_zarr(n, input_axes, t, path)
Construct a Zarr dataset for the level n of a pyramid for the dimensions `input_axes`.
It sets the type to `t` and saves it to `path/n`
"""
function output_zarr(n, input_axes, t, path)
    aggdims = agg_axis.(input_axes, 2^n)
    s = length.(aggdims)
    z = Zeros(t, s...)
    yax = YAXArray(aggdims, z)
    chunked = setchunks(yax , (1024, 1024))
    # This assumes that there is only the spatial dimensions to save
    ds = to_dataset(chunked, )
    dssaved = savedataset(ds; path, skeleton=true, driver=:zarr)
    zar = dssaved.layer.data
    zar
end 

"""
    getpyramids(ras)
Compute the data of the pyramids of a given data cube `ras`.
This returns the data of the pyramids and the dimension values of the aggregated axes.
"""
function getpyramids(reducefunc, ras;recursive=true)
    input_axes = DD.dims(ras)
    n_level = compute_nlevels(ras)            
    if iszero(n_level)
        @info "Array is smaller than the tilesize no pyramids are computed"
        [ras], [dims(ras)]
    end 
    pyramid_sizes =  [ceil.(Int, size(ras) ./ 2^i) for i in 1:n_level]
    pyramid_axes = [agg_axis.(input_axes,2^i) for i in 1:n_level]

    outmin = output_arrays(pyramid_sizes, Float32)
    fill_pyramids(ras,outmin,reducefunc,recursive; threaded=true)

    outmin, pyramid_axes
end

"""
    selectlevel(pyramids, ext, resolution=10; target_imsize=(1024,512)
Internal function to select the raster data that should be plotted on screen. 
`pyramids` is a Vector of Raster objects with increasing coarsity. 
`ext` is the extent of the zoomed in area and `resolution` is the resolution of the data at highest resolution.
`target_imsize` is the target size of the output data that should be plotted.
"""
function selectlevel(pyramid, ext, resolution;target_imsize=(1024, 512))
    # TODO automatically set the target_imsize
    imsize = map(keys(ext)) do bb

        (ext[bb][2] - ext[bb][1]) / resolution[bb]
    end
    dimlevels = log2.(imsize ./ target_imsize)
    minlevel = maximum(dimlevels)
    n_agg = min(max(ceil(Int,minlevel) + 1,1),nlevels(pyramid))
    levels(pyramid)[n_agg][ext]
end


"""
    plot(pyramids)
Plot a Pyramid. 
This will plot the coarsest resolution level at the beginning and will plot higher resolutions after zooming into the plot.
This is expected to be used with interactive Makie backends.
"""
function plot(pyramid::Pyramid;colorbar=true, kwargs...)
    #This should be converted into a proper recipe for Makie but this would depend on a pyramid type.
    fig = Figure()
    lon, lat = DD.dims(parent(pyramid))
    ax = Axis(fig[1,1], limits=(extrema(lon), extrema(lat)), aspect=DataAspect())
    hmap = plot!(ax, pyramid;kwargs...)
    if colorbar
        Colorbar(fig[1,2], hmap)
    end
    ax.autolimitaspect = 1
    FigureAxisPlot(fig, ax, hmap)
end

xkey(keyext) = DD.dim2key(DD.dims(keyext, XDim))
ykey(keyext) = DD.dim2key(DD.dims(keyext, YDim))
#TODO write test, move to utils.jl
function switchkeys(dataext, keyext)
    xk = xkey(keyext)
    yk = ykey(keyext)
    nt = NamedTuple{(xk, yk)}((dataext.X, dataext.Y))
    Extent(nt)
end

function plot!(ax, pyramid::Pyramid;kwargs...)#; rastercrs=crs(parent(pyramid)),plotcrs=EPSG(3857), kwargs...)
    tip = levels(pyramid)[end][:,:]
    #@show typeof(tip)
    data = Observable{DD.AbstractDimMatrix}(tip)
    xval = only(values(Extents.extent(pyramid, XDim)))
    yval = only(values(Extents.extent(pyramid, YDim)))
    rasdataext = Extent(X=xval, Y=yval)
    rasext = extent(pyramid)
    xk = xkey(rasext)
    yk = ykey(rasext)
    resolution = NamedTuple{(xk, yk)}((abs(step(DD.dims(pyramid, XDim))), abs(step(DD.dims(pyramid, YDim)))))
    on(ax.finallimits) do limits
        limext = Extents.extent(limits)
        # Compute limit in raster projection
        #trans = Proj.Transformation(plotcrs, rastercrs, always_xy=true)
        #datalimit = trans_bounds(trans, limext)
        datalimit = switchkeys(limext, rasext)
        if Extents.intersects(rasdataext, limext)
            rasdata = selectlevel(pyramid, datalimit, resolution)
            # Project selected data to plotcrs
            #data.val = Rasters.resample(rasdata, crs=plotcrs, method=:bilinear )
            data.val = rasdata
        end
        notify(data)
    end
    #@show typeof(data)
    hmap = heatmap!(ax, data; interpolate=false, kwargs...)
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

function write(path, pyramid::Pyramid; kwargs...)
    savecube(parent(pyramid), path; kwargs...)
    
    for (i,l) in enumerate(reverse(pyramid.levels))
        outpath = joinpath(path, string(i-1))
        @show outpath
        savecube(l,outpath)
    end
end

"""
    tms_json(dimarray)
Construct a Tile Matrix Set json description from an AbstractDimArray.
This assumes, that we use an ag3ycgregation of two by two pixels in the spatial domain to derive the underlying layers of the pyramids. 
This returns a string representation of the json and is mainly used for writing the TMS definition to the metadata of the Zarr dataset.
"""
function tms_json(pyramid)
    tms = Dict{String, Any}()
    push!(tms, "id"=>"TMS_$(string(DD.name(pyramid)))")
    push!(tms, "title" => "TMS of $(crs(pyramid)) for pyramid")
    push!(tms, "crs" => string(crs(pyramid)))
    push!(tms, "orderedAxes" => pyramidaxes())
    return tms
end

include("broadcast.jl")


end