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
    @test length(axis.scene.plots) == 1
    plot!(axis, pyramid)
    @test length(axis.scene.plots) == 2
end

#= building an RGB pyramid doesn't work, need to think more about it.
@testitem "Pyramid building RGB eltype" begin
    using PyramidScheme: PyramidScheme as PS
    using DimensionalData
    using Colors
    data = rand(RGB, 2000,2000)
    dd = DimArray(data, (X(1:2000), Y(1:2000)))
    pyramid = PS.Pyramid(dd)
    @test pyramid isa PS.Pyramid
    @test eltype(pyramid) isa RGB
end
=#

@testitem "Helper functions" begin
    using PyramidScheme: PyramidScheme as PS
    @test PS.compute_nlevels(zeros(1000)) == 2
    @test PS.compute_nlevels(zeros(1000,1025)) == 3
    @test PS.compute_nlevels(zeros(10000, 8000)) == 6
end

@testitem "ArchGDAL Loading of Geotiff Overviews" begin
    using ArchGDAL: ArchGDAL as AG
    using PyramidScheme: PyramidScheme as PS
    using Rasters
    data = rand(2000,2000)
    r = Raster(data, (X(1:2000), Y(1:2000)))
    tname = tempname() * ".tif"
    write(tname, r, driver="cog", force=true)
    ptif = Pyramid(tname)
    #ras = Raster(path, lazy=true)
    @test ptif isa PS.Pyramid
    @test PS.nlevels(ptif) == 2
    sub = ptif[1:10,1:10]
    @test sub isa PS.Pyramid
end

@testitem "Zarr build Pyramid inplace" begin
    using Zarr
    using PyramidScheme
    using YAXArrays
    using Statistics
    using DimensionalData
    a = rand(1200,1200)
    yax = YAXArray((X(1.:size(a,1)),Y(1.:size(a,2))), a)
    path = tempname() *".zarr"
    savecube(yax, path)
    pyr = buildpyramids(path, resampling_method=mean)
    @test pyr isa Pyramid
    pyrdisk = Pyramid(path)
    @test pyrdisk isa Pyramid 
    @test pyr == pyrdisk
    #pyrmem = PS.Pyramid(yax)
    #@test pyrmem.levels[end][1,1] == pyr.levels[end][1,1]
end

@testitem "selectlevel" begin
    using PyramidScheme: PyramidScheme as PS
    using YAXArrays
    using Extents
    using DimensionalData: DimensionalData as DD
    using DimensionalData.Dimensions
    using Test
    a = rand(1500, 1524)
    yax = YAXArray((X(1.:size(a,1)),Y(1.:size(a,2))), a)
    pyramid = PS.Pyramid(yax)

    target_imsize=(1024, 1024)
    sub =  PS.selectlevel(pyramid, extent(pyramid); target_imsize)
    @test any( (target_imsize ./ 2) .<= size(sub) .<= target_imsize)
    @test size(sub) == (750, 762)
    target_imsize=(256, 256)
    sub =  PS.selectlevel(pyramid, extent(pyramid); target_imsize)
    @test any( (target_imsize ./ 2) .<= size(sub) .<= target_imsize)
    target_imsize=(400, 300)
    sub =  PS.selectlevel(pyramid, extent(pyramid); target_imsize)
    @test any( (target_imsize ./ 2) .<= size(sub) .<= target_imsize)

end


@testitem "cat Pyramids" begin
    using PyramidScheme: PyramidScheme as PS
    using DimensionalData
    pyr1 = Pyramid(rand(X(1:128),Y(1:128)), tilesize=16)
    pyrcat = cat(pyr1, pyr1, dims=Dim{:new}([1,2]))
    pyrcat3 = cat(pyr1, pyr1, pyr1, dims=Dim{:new}([1,2,3]))
    pyr2 = Pyramid(rand(X(129:256), Y(1:128)), tilesize=16)
    pyrcat2 = cat(pyr1, pyr2, dims=X)
    for l in 1:3
        @test PS.levels(pyrcat,l) == cat(PS.levels(pyr1, l), PS.levels(pyr1, l), dims=Dim{:new}([1,2]))
        @test PS.levels(pyrcat3, l) == cat(PS.levels(pyr1, l), PS.levels(pyr1, l), PS.levels(pyr1, l), dims=Dim{:new}([1,2,3]))
        @test PS.levels(pyrcat2, l) == cat(PS.levels(pyr1, l), PS.levels(pyr2, l), dims=X)
    end
end

@testitem "Building pyramid with additional dimension" begin
    # The aim of this test is to check whether we can build a pyramid from a data cube with an extra dimension.
    # We will only build the pyramids on the spatial dimensions and keep the other dimensions as is.
    using YAXArrays
    using Zarr
    using PyramidScheme
    using DimensionalData
    s = (2048, 1024,3)
    a = ones(s)
    yax = YAXArray((X(1.:size(a,1)),Y(1.:size(a,2)), Z(1:3)), a)
    path = tempname() *".zarr"
    savecube(yax, path)
    pyr = buildpyramids(path, resampling_method=sum)
    pyrdisk = Pyramid(path)
    for p in [pyr, pyrdisk]
        @test p isa Pyramid
        @test length(dims(p)) == 3
        @test size(p.levels[end]) == (256,128,3)
        @test p.levels[1][1,1,1] == 4
    end

end
@testitem "Map of pyramids" begin
    using DimensionalData
    pyr1 = Pyramid(fill(1,X(1:128),Y(1:128)), tilesize=16)
    pyr2 = Pyramid(fill(1, X(1:128),Y(1:128)), tilesize=16, resampling_method=sum)
    pyr1_neg = map(x-> x-1, pyr1)
    @test all(all.(iszero, pyr1_neg.levels))
    @test iszero(pyr1_neg.base)
    pyr2_neg = map(x-> x-1, pyr2)
    @test pyr2_neg.levels[1][1,1] == 3

    pyrsum = map((x,y) -> x + y, pyr1, pyr2)
    @test pyrsum[100,30] == 2
    @test pyrsum.levels[1][10,10] == 5
    @test pyrsum.levels[2][10,10] == 17
end

@testitem "Plot of Pyramid" begin
    using CairoMakie
    a = zero(1500, 1524)
    a[1:100, 1:100] .= 2
    yax = YAXArray((X(1.:size(a,1)),Y(1.:size(a,2))), a)
    pyramid = Pyramid(yax)

    # test that colorbar updates
    fig, ax, plt = plot(pyramid)
    cb = Colorbar(fig[1,2], plt.__pyramid_heatmap[])
    @test cb.limits[] == [0.0, 2.0]
    limits!(ax, 200,300, 200, 300)
    @test cb.limits[] == [-0.5, 0.5]
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

