using PyramidScheme:PyramidScheme as PS
using Test
using TestItemRunner
using DimensionalData
using CairoMakie: plot

@run_package_tests
@testitem "Aqua unbound args" begin
    using Aqua
    #Aqua.test_ambiguities([PyramidScheme, Base, Core])
    Aqua.test_unbound_args(PyramidScheme)
end
@testitem "Aqua undefined exports" begin
    using Aqua
    Aqua.test_undefined_exports(PyramidScheme)
end
@testitem "Aqua project extras" begin
        using Aqua
    Aqua.test_project_extras(PyramidScheme)
end

@testitem "Aqua stale deps" begin
        using Aqua
    Aqua.test_stale_deps(PyramidScheme)
end
@testitem "Aqua deps compat" begin
        using Aqua
    Aqua.test_deps_compat(PyramidScheme)
end

@te
@testitem "Pyramid" begin
    using DimensionalData
    using PyramidScheme: PyramidScheme as PS
    using CairoMakie
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
@testitem "Helper functions" begin
    using PyramidScheme: PyramidScheme as PS
    @test PS.compute_nlevels(zeros(1000)) == 2
    @test PS.compute_nlevels(zeros(1000,1025)) == 3
    @test PS.compute_nlevels(zeros(10000, 8000)) == 6
end

@testitem "ArchGDAL Loading of Geotiff Overviews" begin
    using ArchGDAL: ArchGDAL as AG
    using PyramidScheme: PyramidScheme as PS
    @show pwd()
    @show @__DIR__
    path = joinpath(@__DIR__,"data/pyramidmiddle.tif")
    #ras = Raster(path, lazy=true)
    pyr =PS.Pyramid(path)
    @test pyr isa PS.Pyramid
    @test PS.nlevels(pyr) == 2
    sub = pyr[1:10,1:10]
    @test sub isa PS.Pyramid
end

@testitem "Zarr build Pyramid inplace" begin
    using Zarr
    using PyramidScheme: PyramidScheme as PS
    using YAXArrays
    using Statistics
    using DimensionalData
    a = rand(1200,1200)
    yax = YAXArray((X(1.:size(a,1)),Y(1.:size(a,2))), a)
    path = tempname() *".zarr"
    savecube(yax, path)
    pyr = PS.buildpyramids(path, resampling_method=mean)
    @test pyr isa Pyramid
    pyrdisk = Pyramid(path)
    @test pyrdisk isa Pyramid 
    @test pyr .== pyrdisk
    #pyrmem = PS.Pyramid(yax)
    #@test pyrmem.levels[end][1,1] == pyr.levels[end][1,1]
end

#=
@testitem "Comparing zarr pyramid with tif pyramid" begin
    using PyramidScheme: PyramidScheme as PS
    using ArchGDAL
    using Zarr
    using YAXArrays

    pyrtif = PS.Pyramid("test/data/bremen_sea_ice_conc_2022_9_9.tif")
    path = tempname() * ".zarr"
    savecube(parent(pyrtif), path)
    @time pyrzarr = PS.buildpyramids(path)
    pyrdisk = PS.Pyramid(path)
    pyrzarr == pyrtif
end
=#

