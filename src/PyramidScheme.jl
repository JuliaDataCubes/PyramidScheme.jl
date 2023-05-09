module PyramidScheme
using DiskArrayEngine
using DiskArrays
using Zarr
using NetCDF
using Plots
using PyramidScheme

using Statistics
function aggregate_by_factor(xout,x,f)
    fac = ceil(Int,size(x,1)/size(xout,1))
    for j in 1:size(xout,2)
        for i in 1:size(xout,1)
            xout[i,j] = f(view(x,((i-1)*fac+1):min(size(x,1),(i*fac)),((j-1)*fac+1):min(size(x,2),(j*fac))))
        end
    end
end



function all_pyramids!(xout,x,recursive,f)
    xnow = x
    for i in 1:length(xout)
        aggregate_by_factor(xout[i],xnow,f)
        if recursive
            xnow = xout[i]
        end
    end
end
function gen_pyr(x...;recursive=true) 
    allx = Base.front(x)
    f = last(x)
    all_pyramids!(Base.front(allx),last(allx),recursive,f)
end

using DiskArrayEngine: MovingWindow, RegularWindows, InputArray

function fill_pyramids(data, outputs,func,recursive)

    n_level = length(outputs)
    pixel_base_size = 2^n_level
    pyramid_sizes = size.(outputs)
    tmp_sizes = [ceil(Int,pixel_base_size / 2^i) for i in 1:n_level]

    ia  = InputArray(data, windows = arraywindows(size(data),pixel_base_size))

    oa = ntuple(i->create_outwindows(pyramid_sizes[i],windows = arraywindows(pyramid_sizes[i],tmp_sizes[i])),n_level)

    func = create_userfunction(gen_pyr,ntuple(_->eltype(data),length(outputs));is_mutating=true,kwargs = (;recursive),args = (func,))

    op = GMDWop((ia,), oa, func)

    lr = DiskArrayEngine.optimize_loopranges(op,5e8,tol_low=0.2,tol_high=0.05,max_order=2)
    r = DiskArrayEngine.LocalRunner(op,lr,output_arrays,threaded=true)
    run(r)
end

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


function arraywindows(s,w)
    map(s) do l
        RegularWindows(1,l,window=w)
    end
end

using Colors

#testdata = Float32[sin(x)*cos(y) for x in range(0,15,3000), y in range(0,10,2000)];
allatts = NetCDF.open("../../data/C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1.nc","lccs_class").atts

scol = split(allatts["flag_colors"]," ")
colm = [colorant"black";parse.(Color,scol)]

flv = allatts["flag_values"]


testdata = view(NetCDF.open("../../data/C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1.nc","lccs_class"),:,:,1)


lon= ncread("../../data/C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1.nc","lon")
lat= ncread("../../data/C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1.nc","lat")

input_axes = (lon,lat)

n_level = ceil(Int,log2(maximum(size(testdata))/1024))
agg_axis(x,n) = mean.(Iterators.partition(x,n))

pyramid_sizes =  [ceil.(Int, size(testdata) ./ 2^i) for i in 1:n_level]
pyramid_axes = [agg_axis.(input_axes,2^i) for i in 1:n_level]

t = eltype(testdata)
    
function gen_output(t,s)
    outsize = sizeof(t)*prod(s)
    if outsize > 100e6
        p = tempname()
        zcreate(t,s...,path=p,chunks = (1024,1024),fill_value=zero(t))
    else
        zeros(t,s...)
    end
end
output_arrays = [gen_output(UInt8,p) for p in pyramid_sizes]

f = ESALCMode()

fill_pyramids(testdata,output_arrays,f,true)


#Now some plotting code
colmap = similar(colm,256)
for i in 1:length(flv)
    colmap[Int(flv[i])+1] = colm[i]
end
size(testdata)

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



end
