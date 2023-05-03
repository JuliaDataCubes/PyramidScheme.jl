module PyramidScheme
using DiskArrayEngine
using DiskArrays
using Zarr

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


function gen_output(t,s)
    outsize = sizeof(t)*prod(s)
    if outsize > 100e6
        p = tempname()
        zcreate(t,s...,path=p,chunks = (1024,1024),fill_value=zero(t))
    else
        zeros(t,s...)
    end
end


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

end
