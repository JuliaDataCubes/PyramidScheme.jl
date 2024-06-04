module PyramidScheme

using DiskArrayEngine: DiskArrayEngine, MovingWindow, RegularWindows, InputArray, create_outwindows, GMDWop
using DiskArrays: DiskArrays
using Zarr: zcreate
using DimensionalData: DimensionalData as DD
using Extents
using Proj
using Makie: Axis, Colorbar, DataAspect, Figure, FigureAxisPlot, Observable, on, heatmap!
import MakieCore: plot, plot!

using Statistics



"""
    Pyramid
A Pyramid will act as a DimArray on the highest resolution, but subsetting will return another Pyramid.
"""
struct Pyramid{T,N,D,A,B<:DD.AbstractDimArray{T,N,D,A},L} <: DD.AbstractDimArray{T,N,D,A}
    base::B
    levels::L
end

function Pyramid(data::DD.AbstractDimArray)
    pyrdata, pyraxs = getpyramids(mean ∘ skipmissing, data, recursive=false)
    levels = DD.DimArray.(pyrdata, pyraxs)
    Pyramid(data, levels)
end

# refdims
# name

"""
    levels(pyramid::Pyramid)
Return all levels of the `pyramid`. These are order from the base to the coarsest aggregation.
"""
levels(pyramid::Pyramid) = [pyramid.base, pyramid.levels...]
"""
    nlevels(pyramid)
Return the number of levels of the `pyramid`
"""
nlevels(pyramid::Pyramid) = length(levels(pyramid))
Base.parent(pyramid::Pyramid) = pyramid.base
Base.size(pyramid::Pyramid) = size(parent(pyramid))

DD.name(pyramid::Pyramid) = DD.name(parent(pyramid))
DD.refdims(pyramid::Pyramid) = DD.refdims(parent(pyramid))
DD.dims(pyramid::Pyramid) = DD.dims(parent(pyramid))

@inline function DD.rebuild(A::Pyramid, data, dims::Tuple=dims(A), refdims=refdims(A), name=name(A))
    Pyramid(DD.rebuild(parent(A), data, dims, refdims, name, nothing), A.levels)
end

@inline DD.rebuild(A::Pyramid; kw...) = Pyramid(DD.rebuild(parent(A); kw...), A.levels)


@inline function DD.rebuildsliced(f::Function, A::Pyramid, data::AbstractArray, I::Tuple, name=DD.name(A))
    newbase = DD.rebuild(parent(A), parent(data), DD.slicedims(f, A, I)..., name)
    newlevels = map(enumerate(A.levels)) do (z, level)
        Ilevel =  levelindex(z, I)
        leveldata = f(parent(level), Ilevel...)
        DD.rebuild(level, leveldata, DD.slicedims(f, level, Ilevel)..., name)
        
    end
    Pyramid(newbase, newlevels)
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
function fill_pyramids(data, outputs,func,recursive;kwargs...)
    n_level = length(outputs)
    pixel_base_size = 2^n_level
    pyramid_sizes = size.(outputs)
    tmp_sizes = [ceil(Int,pixel_base_size / 2^i) for i in 1:n_level]

    ia  = InputArray(data, windows = arraywindows(size(data),pixel_base_size))

    oa = ntuple(i->create_outwindows(pyramid_sizes[i],windows = arraywindows(pyramid_sizes[i],tmp_sizes[i])),n_level)

    func = DiskArrayEngine.create_userfunction(gen_pyr,ntuple(_->eltype(first(outputs)),length(outputs));is_mutating=true,kwargs = (;recursive),args = (func,))

    op = GMDWop((ia,), oa, func)

    lr = DiskArrayEngine.optimize_loopranges(op,5e8,tol_low=0.2,tol_high=0.05,max_order=2)
    r = DiskArrayEngine.LocalRunner(op,lr,outputs;kwargs...)
    run(r)
end

"""
    ESALCMode(counts)

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
compute_nlevels(data, tilesize=1024) = max(0,ceil(Int,log2(maximum(size(data))/tilesize)))

agg_axis(d,n) = DD.rebuild(d, LinRange(first(d), last(d), cld(length(d), n)))


"""
    gen_output(t,s)

Create output array of type `t` and size `s`
If the array is smaller than 100e6 it is created on disk and otherwise as a temporary zarr file.
"""
function gen_output(t,s)
    outsize = sizeof(t)*prod(s)
    if outsize > 100e6
        p = tempname()
        zcreate(t,s...,path=p,chunks = (1024,1024),fill_value=zero(t))
    else
        zeros(t,s...)
    end
end

"""
    output_arrays(pyramid_sizes)

Create the output arrays for the given `pyramid_sizes`
"""
output_arrays(pyramid_sizes, T) = [gen_output(T,p) for p in pyramid_sizes]


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
function selectlevel(pyramid, ext;target_imsize=(1024, 512))
    # TODO automatically set the target_imsize
    imsize = map(keys(ext)) do bb
        resolution = abs(step(DD.dims(pyramid, bb)))

        (ext[bb][2] - ext[bb][1]) / resolution
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



function plot!(ax, pyramid::Pyramid;kwargs...)#; rastercrs=crs(parent(pyramid)),plotcrs=EPSG(3857), kwargs...)
    data = Observable{DD.AbstractDimArray}(levels(pyramid)[end])
    rasext = extent(pyramid)
    on(ax.finallimits) do limits
        limext = Extents.extent(limits)
        # Compute limit in raster projection
        #trans = Proj.Transformation(plotcrs, rastercrs, always_xy=true)
        #datalimit = trans_bounds(trans, limext)
        datalimit = limext
        if Extents.intersects(rasext, datalimit)
            rasdata = selectlevel(pyramid, datalimit)
            #@show rasdata
            # Project selected data to plotcrs
            #data.val = Rasters.resample(rasdata, crs=plotcrs, method=:bilinear )
            data.val = rasdata
        end
        notify(data)
    end
    hmap = plot!(ax, data; interpolate=false, kwargs...)
end

function trans_bounds(
    trans::Proj.Transformation,
    bbox::Extent,
    densify_pts::Integer = 21,
)::Extent
    xlims, ylims = Proj.bounds(trans, bbox.X, bbox.Y; densify_pts)
    return Extent(X = xlims, Y = ylims)
end

end