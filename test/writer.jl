@testitem "Multiscales store layout" setup = [ConventionHelpers] begin
    using PyramidScheme: PyramidScheme as PS
    using YAXArrays
    using DimensionalData
    using JSON
    using Statistics

    src = YAXArray((X(0.0:49.0), Y(29.0:-1.0:0.0)), rand(Float32, 50, 30))
    dest = tempname() * ".zarr"
    factors = [1, 2, 4]
    pyr = buildpyramids(src, dest; factors, crs="EPSG:4326")
    @test pyr isa Pyramid

    root = JSON.parsefile(joinpath(dest, "zarr.json"))
    @test root["zarr_format"] == 3
    @test root["node_type"] == "group"
    @test !isfile(joinpath(dest, ".zmetadata"))
    expected = PS.multiscales_group_attrs([PS.level_dims(dims(src), f) for f in factors], factors;
        crs="EPSG:4326", resampling_method=mean)
    @test root["attributes"] == group_node(expected)["attributes"]
    s = schemas()
    for conv in ("multiscales", "proj", "spatial")
        @test schema_errors(s[conv], root) === nothing
    end

    # The root only holds metadata, every layout entry is a child group and every child group is in the layout
    layout = root["attributes"]["multiscales"]["layout"]
    assets = [l["asset"] for l in layout]
    @test assets == ["0", "1", "2"]
    @test sort(filter(f -> isdir(joinpath(dest, f)), readdir(dest))) == assets
    consolidated = root["consolidated_metadata"]["metadata"]
    @test sort([k for k in keys(consolidated) if !contains(k, "/")]) == assets
    @test all(consolidated[a]["node_type"] == "group" for a in assets)

    for (level, entry) in zip(assets, layout)
        group = JSON.parsefile(joinpath(dest, level, "zarr.json"))
        @test group["zarr_format"] == 3
        @test group["node_type"] == "group"
        for name in ("layer", "X", "Y")
            array = JSON.parsefile(joinpath(dest, level, name, "zarr.json"))
            @test array["zarr_format"] == 3
            @test array["node_type"] == "array"
            @test haskey(array, "dimension_names")
        end
        data = JSON.parsefile(joinpath(dest, level, "layer", "zarr.json"))
        @test data["dimension_names"] == ["Y", "X"]
        @test data["shape"] == entry["spatial:shape"]
        @test "zstd" in [c["name"] for c in data["codecs"]]
        for name in ("X", "Y")
            coord = JSON.parsefile(joinpath(dest, level, name, "zarr.json"))
            @test "zstd" in [c["name"] for c in coord["codecs"]]
        end
    end

    @test parent(pyr)[:, :] == parent(src)
    pyrmem = Pyramid(src; resampling_method=mean, factors)
    for l in 1:2
        @test Array(parent(PS.levels(pyr, l))) ≈ parent(PS.levels(pyrmem, l))
        @test collect(dims(PS.levels(pyr, l), Y)) ≈ collect(dims(PS.levels(pyrmem, l), Y))
    end
end

@testitem "Dropped edges in cell center style match topozarr" begin
    using PyramidScheme: PyramidScheme as PS
    using YAXArrays
    using DimensionalData
    using JSON

    # Same grid as the odd_size case of fixtures/topozarr/generate.py
    src = YAXArray((Dim{:x}(0.0:49.0), Dim{:y}(29.0:-1.0:0.0)), rand(Float32, 50, 30))
    dest = tempname() * ".zarr"
    pyr = buildpyramids(src, dest; factors=[1, 2, 4], edges=:drop, style=CellCenterStyle(),
        crs="EPSG:4326", spatial_dims=(Dim{:x}, Dim{:y}))
    expected = JSON.parsefile(joinpath(@__DIR__, "fixtures", "topozarr", "odd_size.json"))
    @test JSON.parsefile(joinpath(dest, "zarr.json"))["attributes"] == expected
    @test [size(parent(pyr)); size.(pyr.levels)] == [(50, 30), (25, 15), (12, 7)]
    @test JSON.parsefile(joinpath(dest, "2", "x", "zarr.json"))["shape"] == [12]
    @test JSON.parsefile(joinpath(dest, "2", "y", "zarr.json"))["shape"] == [7]
    @test Array(parent(PS.levels(pyr, 2))) ≈ parent(Pyramid(src; factors=[1, 2, 4], edges=:drop,
        spatial_dims=(Dim{:x}, Dim{:y})).levels[2])
end

@testitem "The source is not modified" begin
    using YAXArrays
    using DimensionalData
    using Zarr

    src = tempname() * ".zarr"
    savecube(YAXArray((X(0.0:63.0), Y(0.0:31.0)), rand(64, 32)), src)
    snapshot(dir) = Dict(relpath(joinpath(root, f), dir) => read(joinpath(root, f))
                         for (root, _, files) in walkdir(dir) for f in files)
    before = snapshot(src)
    dest = tempname() * ".zarr"
    buildpyramids(src, dest; nlevels=2)
    @test snapshot(src) == before

    @test_throws ArgumentError buildpyramids(src, dest; nlevels=2)
    @test buildpyramids(src, dest; nlevels=1, overwrite=true) isa Pyramid
    @test sort(filter(f -> isdir(joinpath(dest, f)), readdir(dest))) == ["0", "1"]
end

@testitem "CRS of the source" begin
    using JSON
    using Rasters
    using DimensionalData

    ras = Raster(rand(Float32, 20, 10), (X(0.0:19.0), Y(9.0:-1.0:0.0)); crs=EPSG(32633))
    dest = tempname() * ".zarr"
    buildpyramids(ras, dest; factors=[1, 2])
    attrs = JSON.parsefile(joinpath(dest, "zarr.json"))["attributes"]
    @test attrs["proj:code"] == "EPSG:32633"
    @test "proj" in [c["name"] for c in attrs["zarr_conventions"]]

    buildpyramids(ras, dest; factors=[1, 2], crs=nothing, overwrite=true)
    attrs = JSON.parsefile(joinpath(dest, "zarr.json"))["attributes"]
    @test !any(startswith("proj:"), keys(attrs))
    @test [c["name"] for c in attrs["zarr_conventions"]] == ["multiscales", "spatial"]
end

@testitem "Multiscales store with an extra dimension" setup = [ConventionHelpers] begin
    using PyramidScheme: PyramidScheme as PS
    using YAXArrays
    using DimensionalData
    using JSON

    src = YAXArray((X(0.0:31.0), Y(0.0:15.0), Ti(1:3)), rand(32, 16, 3))
    dest = tempname() * ".zarr"
    pyr = buildpyramids(src, dest; factors=[1, 2], chunksize=8)
    root = JSON.parsefile(joinpath(dest, "zarr.json"))
    @test root["attributes"]["multiscales"]["layout"][2]["transform"]["scale"] == [1.0, 2.0, 2.0]
    data = JSON.parsefile(joinpath(dest, "1", "layer", "zarr.json"))
    @test data["dimension_names"] == ["Ti", "Y", "X"]
    @test data["shape"] == [3, 8, 16]
    @test data["chunk_grid"]["configuration"]["chunk_shape"] == [1, 8, 8]
    @test schema_errors(schemas()["multiscales"], root) === nothing
    @test size(pyr.levels[1]) == (16, 8, 3)
end

@testitem "Write an in memory pyramid" setup = [ConventionHelpers] begin
    using PyramidScheme: PyramidScheme as PS
    using YAXArrays
    using DimensionalData
    using JSON

    src = DimArray(rand(40, 24), (X(0.0:39.0), Y(23.0:-1.0:0.0)))
    pyrmem = Pyramid(src; resampling_method=sum, factors=[1, 2, 4])
    dest = tempname() * ".zarr"
    pyr = write(dest, pyrmem; crs="EPSG:4326")
    root = JSON.parsefile(joinpath(dest, "zarr.json"))
    @test root["attributes"]["multiscales"]["resampling_method"] == "sum"
    @test [l["asset"] for l in root["attributes"]["multiscales"]["layout"]] == ["0", "1", "2"]
    s = schemas()
    for conv in ("multiscales", "proj", "spatial")
        @test schema_errors(s[conv], root) === nothing
    end
    for l in 0:2
        @test Array(parent(PS.levels(pyr, l))) ≈ parent(PS.levels(pyrmem, l))
    end
    @test_throws ArgumentError write(tempname() * ".zarr", Pyramid(src, pyrmem.levels, Dict()))
end
