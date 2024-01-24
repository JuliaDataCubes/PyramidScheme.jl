module PyramidScheme

using DiskArrayEngine: DiskArrayEngine, MovingWindow, RegularWindows, InputArray, create_outwindows, GMDWop
using DiskArrays: DiskArrays
using Zarr: zcreate
using DimensionalData: rebuild, dims
using GLMakie
using Rasters
using Extents
using Proj

using Statistics

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

agg_axis(x,n) = rebuild(x, mean.(Iterators.partition(x,n)))


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
    input_axes = (dims(ras))
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
function selectlevel(pyramids, ext, resolution=10;target_imsize=(1024, 512))
    imsize = (ext.X[2] - ext.X[1], ext.Y[2] - ext.Y[1]) ./ resolution
    n_agg = min(max(ceil(Int,minimum(log2.(imsize ./ target_imsize))) + 1,1),length(pyramids))
    pyramids[n_agg][ext]
end


"""
    plotpyramids(pyramids)
A helper function to plot a dataset of pyramids.
At the moment pyramids is expected to be a list of pyramids with the same extent. 
"""
function plotpyramids(pyramids;colorbar=true, kwargs...)
    #This should be converted into a proper recipe for Makie but this would depend on a pyramid type.
    fig = Figure()
    lon, lat = dims(pyramids[1])
    ax = Axis(fig[1,1], limits=(extrema(lon), extrema(lat)), aspect=DataAspect())
    hmap = plotpyramids!(ax, pyramids;kwargs...)
    if colorbar
        Colorbar(fig[1,2], hmap)
    end
    ax.autolimitaspect = 1
    fig, ax, hmap
end



function plotpyramids!(ax, pyramids;rastercrs=crs(pyramids[1]),plotcrs=EPSG(3857), kwargs...)
    data = Observable{Raster}(pyramids[end])
    rasext = extent(pyramids[end])
    on(ax.finallimits) do limits
        limext = Extents.extent(limits)
        # Compute limit in raster projection
        trans = Proj.Transformation(plotcrs, rastercrs, always_xy=true)
        datalimit = trans_bounds(trans, limext)

        if Extents.intersects(rasext, datalimit)
            rasdata = selectlevel(pyramids, datalimit)
            # Project selected data to plotcrs
            data.val = Rasters.resample(rasdata, crs=plotcrs)
        end
        notify(data)
    end
    
    hmap = heatmap!(ax, data; interpolate=false, kwargs...)
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