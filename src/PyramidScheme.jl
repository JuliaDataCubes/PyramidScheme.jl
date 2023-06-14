module PyramidScheme
using DiskArrayEngine: DiskArrayEngine, MovingWindow, RegularWindows, InputArray, create_outwindows, GMDWop
using DiskArrays: DiskArrays
using Zarr: zcreate
using DimensionalData: rebuild

using Statistics

"""
    aggregate_by_factor(xout, x, f)

Aggregate the data `x` using the function `f` and store the results in `xout`.
This function is used as the inner function for the DiskArrayEngine call.
"""
function aggregate_by_factor(xout,x,f)
    fac = ceil(Int,size(x,1)/size(xout,1))
    for j in 1:size(xout,2)
        for i in 1:size(xout,1)
            xout[i,j] = f(view(x,((i-1)*fac+1):min(size(x,1),(i*fac)),((j-1)*fac+1):min(size(x,2),(j*fac))))
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
    for i in 1:length(xout)
        @debug "Pyramidnumber $i"
        aggregate_by_factor(xout[i],xnow,f)
        if recursive
            xnow = xout[i]
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
compute_nlevels(data, tilesize=1024) = ceil(Int,log2(maximum(size(data))/tilesize))

agg_axis(x,n) = rebuild(x, mean.(Iterators.partition(x,n)))

#pyramid_sizes(n_level) =  [ceil.(Int, size(testdata) ./ 2^i) for i in 1:n_level]
#pyramid_axes(n_level) = [agg_axis.(input_axes,2^i) for i in 1:n_level]

#t = eltype(testdata)

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
    selectlevel(pyramids, ext, resolution=10; target_imsize=(1024,512)
Internal function to select the raster data that should be plotted on screen. 
`pyramids` is a Vector of Raster objects with increasing coarsity. 
`ext` is the extent of the zoomed in area and `resolution` is the resolution of the data at highest resolution.
`target_imsize` is the target size of the output data that should be plotted.
"""
function selectlevel(pyramids, ext, resolution=10;target_imsize=(1024, 512))
    imsize = (ext.X[2] - ext.X[1], ext.Y[2] - ext.Y[1]) ./ resolution
    @show imsize
    n_agg = min(max(ceil(Int,minimum(log2.(imsize ./ target_imsize))) + 1,1),length(pyramids))
    @show n_agg
    pyramids[n_agg][ext]
end
#f = ESALCMode()

#fill_pyramids(testdata,output_arrays,f,true)


#Now some plotting code
#colmap = similar(colm,256)
#for i in 1:length(flv)
#    colmap[Int(flv[i])+1] = colm[i]
#end
#size(testdata)


#=
function plot_im(data,colmap,cent,imsize;target_imsize = (1024,512))
    n_agg = ceil(Int,minimum(log2.(imsize./target_imsize)))+1
    n_agg = clamp(n_agg,1,length(data))
    @show n_agg
    b1 = (cent .- imsize .÷  2 .+1) .÷ (2^(n_agg-1))
    b2 = (cent .+ imsize .÷ 2) .÷ (2^(n_agg-1))
    datanow = data[n_agg]
    snow = size(datanow)
    b1c = clamp.(b1,1,snow)
    b2c = clamp.(b2,1,snow)
    inds_intern = range.(b1c,b2c)
    x = data[n_agg][inds_intern...]
    permutedims(map(x) do ix
        colmap[Int(ix)+1]
    end)
end


mutable struct Viewer
    data
    colmap
    center
    zoom
end
zoom_out(v::Viewer) = v.zoom = v.zoom*1.3
zoom_in(v::Viewer) = v.zoom = v.zoom/1.3
move_right(v::Viewer) = v.center = round.(Int,(v.center[1]+300*v.zoom,v.center[2]))
move_left(v::Viewer) = v.center = round.(Int,(v.center[1]-300*v.zoom,v.center[2]))
move_down(v::Viewer) = v.center = round.(Int,(v.center[1],v.center[2]+300*v.zoom))
move_up(v::Viewer) = v.center = round.(Int,(v.center[1],v.center[2]-300*v.zoom))
function doplot(v::Viewer)
    target_imsize = (1024,512)
    imsize = round.(Int,v.zoom .*target_imsize)
    plot_im(v.data,v.colmap, v.center,imsize)
end

function apply_input(v,x)
    if x=='a'
        move_left(v)
    elseif x == 'd'
        move_right(v)
    elseif x == 'w'
        move_up(v)
    elseif x == 's'
        move_down(v)
    elseif x == 'k'
        zoom_in(v)
    elseif x == 'l'
        zoom_out(v)
    else 
        return 1
    end
    0
end

data = [testdata,output_arrays...]

v = Viewer(data,colmap,(80000,30000),1.0)

while true
    inp = readline(stdin)
    isempty(inp) && break
    r = apply_input(v,first(inp))
    if r == 1
        break
    else
        display(doplot(v))
    end
end

=#

end
