@testitem "Read multiscale stores" begin
    using PyramidScheme: PyramidScheme as PS
    using YAXArrays
    using DimensionalData

    src = YAXArray((X(0.0:49.0), Y(29.0:-1.0:0.0)), rand(Float32, 50, 30))
    for (factors, edges) in (([1, 2, 4], :keep), ([1, 2, 8], :drop))
        dest = tempname() * ".zarr"
        pyr = buildpyramids(src, dest; factors, edges, crs="EPSG:4326")
        disk = Pyramid(dest)
        @test disk isa Pyramid
        @test disk == pyr
        @test PS.level_factors(disk) == factors
        @test metadata(disk)["proj:code"] == "EPSG:4326"
        @test haskey(metadata(disk), "multiscales")
    end
end

@testitem "Read a multiscales layout in any order" begin
    using PyramidScheme: PyramidScheme as PS
    using YAXArrays
    using DimensionalData
    using JSON

    src = YAXArray((X(0.0:31.0), Y(0.0:31.0)), rand(32, 32))
    dest = tempname() * ".zarr"
    pyr = buildpyramids(src, dest; factors=[1, 2, 4], consolidated=false)
    rootfile = joinpath(dest, "zarr.json")
    root = JSON.parsefile(rootfile)
    reverse!(root["attributes"]["multiscales"]["layout"])
    write(rootfile, JSON.json(root))
    @test Pyramid(dest) == pyr
end

@testitem "Scales of a multiscales layout" begin
    using PyramidScheme: PyramidScheme as PS
    using JSON

    # Official example with downsampled and upsampled levels
    example = JSON.parsefile(joinpath(@__DIR__, "schemas", "examples", "dem-multiresolution.json"))
    layout = example["attributes"]["multiscales"]["layout"]
    @test [l["asset"] for l in layout] == ["dem_30m", "dem_90m", "dem_270m", "dem_10m_superres"]
    @test PS.layout_scales(layout) ≈ [1.0, 3.0, 9.0, 0.333]

    # Non spatial dimensions have a scale of one, level 0 may have no transform
    withtime = [Dict("asset" => "0"),
        Dict("asset" => "1", "derived_from" => "0", "transform" => Dict("scale" => [1.0, 2.0, 2.0]))]
    @test PS.layout_scales(withtime) == [1.0, 2.0]
    missing_source = [Dict("asset" => "0"), Dict("asset" => "1", "derived_from" => "2")]
    @test_throws ArgumentError PS.layout_scales(missing_source)
    cycle = [Dict("asset" => "0", "derived_from" => "1", "transform" => Dict("scale" => [2.0])),
        Dict("asset" => "1", "derived_from" => "0", "transform" => Dict("scale" => [2.0]))]
    @test_throws ArgumentError PS.layout_scales(cycle)
end

@testitem "Read pyramids written by older versions" begin
    using PyramidScheme: PyramidScheme as PS
    using YAXArrays
    using DimensionalData
    using Zarr

    # Base at the root and the levels in the groups 1, 2, ... next to it, as zarr v2
    base = YAXArray((X(0.0:31.0), Y(0.0:15.0)), rand(32, 16))
    mem = Pyramid(base; factors=[1, 2, 4])
    path = tempname() * ".zarr"
    savecube(base, path)
    for (i, level) in enumerate(mem.levels)
        savecube(YAXArray(dims(level), parent(level)), joinpath(path, string(i)))
    end
    pyr = Pyramid(path)
    @test PS.nlevels(pyr) == 2
    @test size.(pyr.levels) == [(16, 8), (8, 4)]
    @test Array(parent(pyr.levels[2])) ≈ parent(mem.levels[2])
end

@testitem "Read a store written by topozarr" begin
    using PyramidScheme: PyramidScheme as PS
    using DimensionalData
    using Statistics

    # Sharded arrays, xarray coordinates and a scalar spatial_ref array without dimension_names.
    # topozarr 0.1.9 measured the translation between levels between cell centers, 0.2.0 writes zero.
    for store in ("store.zarr", "store-0.1.9.zarr")
        path = joinpath(@__DIR__, "fixtures", "topozarr", store)
        pyr = Pyramid(path)
        @test pyr isa Pyramid
        @test size(parent(pyr)) == (32, 32)
        @test size.(pyr.levels) == [(16, 16)]
        @test name.(dims(parent(pyr))) == (:x, :y)
        base, level = Array(parent(parent(pyr))), Array(parent(pyr.levels[1]))
        # Values written by generate.py: v[y, x] = 32 * y + x, zero based
        @test base[3, 2] == 32 + 2
        @test all(level[i, j] ≈ mean(base[2i-1:2i, 2j-1:2j]) for i in 1:16, j in 1:16)
        @test collect(lookup(pyr.levels[1], Dim{:x}))[1:2] == [0.5, 2.5]
        @test metadata(pyr)["proj:code"] == "EPSG:4326"
    end
end

@testitem "Select levels with arbitrary factors" begin
    using PyramidScheme: PyramidScheme as PS
    using DimensionalData
    using Extents

    pyr = Pyramid(DimArray(rand(1024, 1024), (X(1.0:1024.0), Y(1.0:1024.0))); factors=[1, 4, 8])
    # 1024 / 400 = 2.56 cells per pixel: the level downsampled by 4 is the finest one that fits
    @test size(PS.selectlevel(pyr, extent(pyr); target_imsize=(400, 400))) == (256, 256)
    @test size(PS.selectlevel(pyr, extent(pyr); target_imsize=(2048, 2048))) == (1024, 1024)
    @test size(PS.selectlevel(pyr, extent(pyr); target_imsize=(50, 50))) == (128, 128)
end
