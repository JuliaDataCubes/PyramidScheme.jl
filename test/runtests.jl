using PyramidScheme:PyramidScheme as PS
using Test
using DimensionalData
#using Zarr
#using ArchGDAL
#using NetCDF
#using YAXArrayBase
#using Rasters
#using Statistics
#using YAXArrays
#using Extents


@testset "Pyramid" begin
    using DimensionalData
    using PyramidScheme: PyramidScheme as PS
    data = rand(2000,2000)
    dd = DimArray(data, (X(1:2000), Y(1:2000)))
    pyramid = PS.Pyramid(dd)
    @test PS.nlevels(pyramid) == 2
    subpyramid = pyramid[X=1..10, Y=1..10]
    @test subpyramid isa PS.Pyramid
    @test size(subpyramid) == (10,10)
    @test parent(subpyramid) == data[1:10,1:10]
end
@testset "Helper functions" begin
    @test PS.compute_nlevels(zeros(1000)) == 0
    @test PS.compute_nlevels(zeros(1000,1025)) == 1
    @test PS.compute_nlevels(zeros(10000, 8000)) == 4
end

@testset "User facing functions" begin
    data = rand(256,256)
    dd = DimArray(data, (X(1:256), Y(1:256)))
    pyramid = Pyramid(dd, tilesize=8)
    @test pyr[1] == dd
    @test pyr isa AbstractDimArray

    pyramid[level=1]
    @test plot(pyr)
    

end


#testdata = Float32[sin(x)*cos(y) for x in range(0,15,3000), y in range(0,10,2000)];
cd("test")
inpath = "/home/fcremer/Documents/rqa_cscale/data/ForestType/2017_FOREST-CLASSES_EU010M_E048N018T1.tif"
inpath = "/home/fcremer/Documents/rqa_cscale/data/ger"
inpath = "data/E048N018T3_rqatrend_VH_A117__thresh_4_8_8_2023-03-22T16:00:28.438.zarr/"


#arr = zopen(inpath,fill_as_missing=true)
#allatts = NetCDF.open("../../data/C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1.nc","lccs_class").atts
c = Cube(inpath)
ras = Raster(yaxconvert(DimArray, c))

#lon = arr.arrays["X"]
#lat = arr.arrays["Y"]
#scol = split(allatts["flag_colors"]," ")
#colm = [colorant"black";parse.(Color,scol)]

#flv = allatts["flag_values"]
#input_axes = (1:100,1:100)

#input_axes = (c.X, c.Y)
input_axes = (dims(ras))
#testdata = arr.arrays["layer"]
n_level = PS.compute_nlevels(c)
pyramid_sizes =  [ceil.(Int, size(c) ./ 2^i) for i in 1:n_level]
pyramid_axes = [PS.agg_axis.(input_axes,2^i) for i in 1:n_level]

outarrs = PS.output_arrays(pyramid_sizes, Float32)
function myfunc(x)
    x2 = filter(i-> i<2e4 && !isnan(i),x)
    if !isempty(x2)
        mean(x2)
    else
        NaN32
    end
end



#@time PS.fill_pyramids(c.data[:,:],outarrs,x->sum(x) >0,true)
@time PS.fill_pyramids(c.data[:,:],outarrs,mean,true)

pyramids = [ras, Raster.(outarrs, pyramid_axes)...]

t = eltype(testdata)

#testdata = view(NetCDF.open("../../data/C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1.nc","lccs_class"),:,:,1)

#=
function plot_im(data,xlim,ylim,pyramid_axes;target_imsize = (1024,512))

    imsize = (xlim[2]-xlim[1]),(ylim[2]-ylim[1])
    @show imsize
    n_agg = ceil(Int,minimum(log2.(imsize./target_imsize))) + 1
    @show n_agg
    n_agg = clamp(n_agg,1,length(data))
    b1 = round.(Int,xlim) .÷ (2^(n_agg-1))
    @show b1
    b2 = round.(Int,ylim) .÷ (2^(n_agg-1))
    datanow = data[n_agg]
    snow = size(datanow)
    b1c = clamp.(b1,1,snow)
    b2c = clamp.(b2,1,snow)
    inds_intern = range(b1c...), range(b2c...)
    @show b1c
    @show inds_intern
    xvals = pyramid_axes[n_agg][1][inds_intern[1]]
    #@show xvals
    yvals = pyramid_axes[n_agg][2][inds_intern[2]]
    arr = data[n_agg][inds_intern...]
    @show size.([arr, xvals, yvals])
    xvals, yvals, xlim,ylim,arr
end
=#
#plot_im(outarrs, xlim, ylim,pyramid_axes)

#lon = c.X
#lat = c.Y
using GLMakie


lon,lat = input_axes
#xvals = Observable(collect(Float64, input_axes[1]))
#yvals = Observable(collect(Float64, input_axes[2]))
#xlim = Observable(Float64.(extrema(lon)))
#ylim = Observable(Float64.(extrema(lat)))
fig = Figure()

ax = Axis(fig[1,1], limits=(extrema(lon), extrema(lat)))
data = Observable{Raster}(pyramids[end])

#image!(ax, ras)
#(pyramid_axes[end]...,outarrs[end])
on(ax.finallimits) do limits
    data.val = PS.selectlevel(pyramids, extent(limits))
    notify(data)
end
@show size(data.val)
#xvals.v = lift(i-> i[1], data)
#yvals = @lift $data[2]
#xlimob = @lift $data[3]
#ylimob = @lift $data[4]
#arr = @lift $data[5]

image!(ax, data, interpolate=false)

histvals =@lift filter(isfinite, skipmissing($data))
axhist = Axis(fig[2,1])
hist!(axhist, histvals)

#heatmap!(ax, xvals, yvals, ras,xlim=xlim,ylim=ylim)
