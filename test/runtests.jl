using PyramidScheme:PyramidScheme as PS
using Test
using DimensionalData
using CairoMakie: plot

@testset "Pyramid" begin
    using DimensionalData
    using PyramidScheme: PyramidScheme as PS
    data = zeros(2000,2000)
    dd = DimArray(data, (X(1:2000), Y(1:2000)))
    pyramid = PS.Pyramid(dd)
    #@test PS.nlevels(pyramid) == 2
    subpyramid = pyramid[X=1..10, Y=1..10]
    @test subpyramid isa PS.Pyramid
    @test size(subpyramid) == (10,10)
    @test parent(subpyramid) == data[1:10,1:10]
    fig, axis, h = plot(pyramid)
end
@testset "Helper functions" begin
    @test PS.compute_nlevels(zeros(1000)) == 0
    @test PS.compute_nlevels(zeros(1000,1025)) == 1
    @test PS.compute_nlevels(zeros(10000, 8000)) == 4
end

@testset "ArchGDAL Loading of Geotiff Overviews" begin
    using ArchGDAL: ArchGDAL as AG
    using PyramidScheme: PyramidScheme as PS
    using Rasters

    path = "test/data/pyramidmiddle.tif"
    #ras = Raster(path, lazy=true)
    pyr =PS.Pyramid(path)
    @test pyr isa PS.Pyramid
    plot(pyr)
    @test PS.nlevels(pyr) == 3
    sub = pyr[1:10,1:10]
    @test sub isa Pyramid
end