@testmodule ConventionHelpers begin
using DimensionalData: DimensionalData as DD
using JSON: JSON
using JSONSchema: Schema, validate
using PyramidScheme: PyramidScheme as PS

export schemas, conventionschema, group_node, schema_errors, coarsen_dim, pyramid_dims

const SCHEMADIR = joinpath(@__DIR__, "schemas")

conventionschema(name) = JSON.parsefile(joinpath(SCHEMADIR, "$name-v0.1.json"))

# The proj schema references the remote PROJJSON schema, which JSONSchema.jl downloads.
# proj:projjson is never written, so the reference is replaced to keep the tests offline.
function schemas()
    proj = conventionschema("proj")
    proj["\$defs"]["projAttributes"]["properties"]["proj:projjson"] = Dict("type" => "object")
    Dict("multiscales" => Schema(conventionschema("multiscales")),
        "proj" => Schema(proj),
        "spatial" => Schema(conventionschema("spatial")))
end

# zarr.json of a group with these attributes, round tripped through JSON as it would be on disk
group_node(attrs) = JSON.parse(JSON.json(Dict("zarr_format" => 3, "node_type" => "group", "attributes" => attrs)))

schema_errors(schema, node) = validate(schema, node)

# Coarsen a point lookup by `factor`, keeping partial windows at the end
function coarsen_dim(d, factor)
    res = PS.cellsize(d)
    corner = first(DD.lookup(d)) - res / 2
    newres = res * factor
    DD.basetypeof(d)(range(corner + newres / 2; step=newres, length=cld(length(d), factor)))
end

pyramid_dims(dims, factors; spatial_dims=(DD.XDim, DD.YDim)) =
    [map(d -> DD.hasdim(DD.dims(dims, spatial_dims), d) ? coarsen_dim(d, f) : d, dims) for f in factors]
end

@testitem "Convention entries match the schemas" setup = [ConventionHelpers] begin
    using PyramidScheme: PyramidScheme as PS
    for (conv, name) in [(PS.MULTISCALES_CONVENTION, "multiscales"), (PS.PROJ_CONVENTION, "proj"),
        (PS.SPATIAL_CONVENTION, "spatial")]
        consts = conventionschema(name)["\$defs"]["conventionMetadata"]["properties"]
        attrs = PS.convention_attrs(conv)
        @test sort(collect(keys(attrs))) == sort(collect(keys(consts)))
        for (k, v) in attrs
            @test v == consts[k]["const"]
        end
    end
end

@testitem "Spatial transform" setup = [ConventionHelpers] begin
    using PyramidScheme: PyramidScheme as PS
    using DimensionalData
    using DimensionalData.Lookups

    up =(X(0.0:15.0), Y(0.0:15.0))
    @test PS.spatial_transform(up...) == [1.0, 0.0, -0.5, 0.0, 1.0, -0.5]
    down = (X(0.0:15.0), Y(15.0:-1.0:0.0))
    @test PS.spatial_transform(down...) == [1.0, 0.0, -0.5, 0.0, -1.0, 15.5]
    @test PS.spatial_attrs(up...)["spatial:bbox"] == PS.spatial_attrs(down...)["spatial:bbox"] == [-0.5, -0.5, 15.5, 15.5]

    attrs = PS.spatial_attrs(X(10.0:2.0:18.0), Y(0.0:3.0:6.0))
    @test attrs["spatial:dimensions"] == ["Y", "X"]
    @test attrs["spatial:shape"] == [3, 5]
    @test attrs["spatial:registration"] == "pixel"
    @test attrs["spatial:bbox"] == [9.0, -1.5, 19.0, 7.5]
    # Names are taken from the dimensions, the order is always [y, x]
    lonlat = (Dim{:lon}(0.0:1.0:3.0), Dim{:lat}(0.0:1.0:1.0))
    @test PS.spatial_attrs(lonlat...)["spatial:dimensions"] == ["lat", "lon"]

    # Interval lookups are converted to cell centers according to their locus
    start = X(Sampled(0.0:10.0:20.0; sampling=Intervals(Start()), order=ForwardOrdered(), span=Regular(10.0)))
    stop = X(Sampled(10.0:10.0:30.0; sampling=Intervals(End()), order=ForwardOrdered(), span=Regular(10.0)))
    center = X(Sampled(5.0:10.0:25.0; sampling=Intervals(Center()), order=ForwardOrdered(), span=Regular(10.0)))
    for d in (start, stop, center)
        @test PS.spatial_transform(d, Y(0.0:1.0))[1:3] == [10.0, 0.0, 0.0]
    end

    @test_throws ArgumentError PS.spatial_transform(X([0.0, 1.0, 3.0]), Y(0.0:1.0))
    @test_throws ArgumentError PS.spatial_transform(X([5.0]), Y(0.0:1.0))
    @test PS.spatial_transform(X([5.0]), Y(0.0:1.0); fallback_res=(4.0, 1.0))[1:3] == [4.0, 0.0, 3.0]
end

@testitem "proj attributes" setup = [ConventionHelpers] begin
    using PyramidScheme: PyramidScheme as PS
    using Rasters: EPSG
    Proj = PS.Proj

    attrs = PS.proj_attrs("EPSG:4326")
    @test attrs["proj:code"] == "EPSG:4326"
    @test startswith(attrs["proj:wkt2"], "GEOGCRS")
    @test PS.proj_attrs(EPSG(32633))["proj:code"] == "EPSG:32633"
    # A WKT2 string that carries an ID keeps its code
    wkt = Proj.proj_as_wkt(Proj.CRS("EPSG:3857"), Proj.PJ_WKT2_2019)
    @test PS.proj_attrs(wkt)["proj:code"] == "EPSG:3857"
    # Without an unambiguous authority code only the WKT2 is written
    noid = PS.proj_attrs("+proj=longlat +datum=WGS84 +no_defs +type=crs")
    @test !haskey(noid, "proj:code")
    @test haskey(noid, "proj:wkt2")
end

@testitem "Resampling method names" begin
    using PyramidScheme: PyramidScheme as PS
    using Statistics
    @test PS.resampling_method_name(mean) == "average"
    @test PS.resampling_method_name(mean ∘ skipmissing) == "average"
    @test PS.resampling_method_name(median) == "med"
    @test PS.resampling_method_name(sum) == "sum"
    @test PS.resampling_method_name(maximum) == "max"
    @test PS.resampling_method_name(minimum) == "min"
    @test PS.resampling_method_name(first) == "nearest"
    @test PS.resampling_method_name(PS.ESALCMode()) == "mode"
    @test PS.resampling_method_name(std) == "std"
end

@testitem "Multiscales layout" setup = [ConventionHelpers] begin
    using PyramidScheme: PyramidScheme as PS
    using DimensionalData
    using Statistics

    dims = (X(0.0:31.0), Y(31.0:-1.0:0.0), Ti(1:3))
    levels = pyramid_dims(dims, [1, 2, 4])
    ms = PS.multiscales_attrs(levels, [1, 2, 4]; resampling_method=mean)
    layout = ms["layout"]
    @test ms["resampling_method"] == "average"
    @test [l["asset"] for l in layout] == ["0", "1", "2"]
    @test !haskey(layout[1], "derived_from")
    @test !haskey(layout[1], "resampling_method")
    @test [l["derived_from"] for l in layout[2:end]] == ["0", "1"]
    @test all(l["resampling_method"] == "average" for l in layout[2:end])
    # Stored order is (Ti, Y, X)
    @test layout[1]["transform"] == Dict("scale" => [1.0, 1.0, 1.0], "translation" => [0.0, 0.0, 0.0])
    @test layout[2]["transform"] == Dict("scale" => [1.0, 2.0, 2.0], "translation" => [0.0, 0.0, 0.0])
    @test layout[3]["transform"]["scale"] == [1.0, 2.0, 2.0]
    @test [l["spatial:shape"] for l in layout] == [[32, 32], [16, 16], [8, 8]]
    t = [l["spatial:transform"] for l in layout]
    @test [ti[1] for ti in t] == [1.0, 2.0, 4.0]
    @test [ti[5] for ti in t] == [-1.0, -2.0, -4.0]
    @test all(ti[[3, 6]] == [-0.5, 31.5] for ti in t)

    # Levels computed from the base are derived from level 0
    nonrec = PS.multiscales_attrs(levels, [1, 2, 4]; resampling_method=mean, recursive=false)["layout"]
    @test [l["derived_from"] for l in nonrec[2:end]] == ["0", "0"]
    @test nonrec[3]["transform"]["scale"] == [1.0, 4.0, 4.0]

    # Sparse levels use the ratio to their source
    sparse = PS.multiscales_attrs(pyramid_dims(dims, [1, 2, 8]), [1, 2, 8]; resampling_method=sum)["layout"]
    @test sparse[3]["transform"]["scale"] == [1.0, 4.0, 4.0]
    @test sparse[3]["resampling_method"] == "sum"

    # Partial windows at the end are kept, the origin stays the same
    odd = pyramid_dims((X(1.0:1500.0), Y(1.0:1524.0)), [1, 2, 4])
    oddlayout = PS.multiscales_attrs(odd, [1, 2, 4]; resampling_method=mean)["layout"]
    @test [l["spatial:shape"] for l in oddlayout] == [[1524, 1500], [762, 750], [381, 375]]
    @test all(l["spatial:transform"][[3, 6]] == [0.5, 0.5] for l in oddlayout)

    # A level collapsed to one cell takes its cell size from level 0
    tiny = pyramid_dims((X(0.0:2.0), Y(0.0:7.0)), [1, 4])
    @test PS.multiscales_attrs(tiny, [1, 4]; resampling_method=mean)["layout"][2]["spatial:transform"][1] == 4.0

    @test_throws ArgumentError PS.multiscales_attrs(levels, [1, 2]; resampling_method=mean)
    for bad in ([2, 4, 8], [1, 4, 2], [1, 2, 3], Int[])
        @test_throws ArgumentError PS.check_factors(bad)
    end
end

@testitem "Multiscales group attributes validate against the schemas" setup = [ConventionHelpers] begin
    using PyramidScheme: PyramidScheme as PS
    using DimensionalData
    using Statistics

    s = schemas()
    cases = [
        ("square", (X(0.0:31.0), Y(0.0:31.0)), [1, 2, 4], "EPSG:4326"),
        ("non square", (X(0.0:47.0), Y(0.0:31.0)), [1, 2, 4], "EPSG:4326"),
        ("north up", (X(0.0:31.0), Y(31.0:-1.0:0.0)), [1, 2, 4], "EPSG:4326"),
        ("extra dimension", (X(0.0:31.0), Y(0.0:31.0), Ti(1:3)), [1, 2, 4], "EPSG:4326"),
        ("projected", (X(0.0:100.0:3100.0), Y(3100.0:-100.0:0.0), Ti(1:2)), [1, 2], "EPSG:3857"),
        ("sparse factors", (X(0.0:63.0), Y(0.0:47.0)), [1, 2, 8], "EPSG:4326"),
        ("odd size", (X(1.0:1500.0), Y(1.0:1524.0)), [1, 2, 4], "EPSG:4326"),
    ]
    for (name, dims, factors, crs) in cases
        attrs = PS.multiscales_group_attrs(pyramid_dims(dims, factors), factors; crs, resampling_method=mean)
        node = group_node(attrs)
        for conv in ("multiscales", "proj", "spatial")
            @test schema_errors(s[conv], node) === nothing
        end
        @test [c["name"] for c in node["attributes"]["zarr_conventions"]] == ["multiscales", "proj", "spatial"]
        @test node["attributes"]["proj:code"] == crs
    end

    # Without a CRS the proj convention is left out
    attrs = PS.multiscales_group_attrs(pyramid_dims((X(0.0:31.0), Y(0.0:31.0)), [1, 2]), [1, 2]; resampling_method=mean)
    node = group_node(attrs)
    @test schema_errors(s["multiscales"], node) === nothing
    @test schema_errors(s["spatial"], node) === nothing
    @test [c["name"] for c in node["attributes"]["zarr_conventions"]] == ["multiscales", "spatial"]
    @test !any(startswith("proj:"), keys(node["attributes"]))

    # The schemas reject what they should
    bad = deepcopy(node)
    delete!(bad["attributes"]["multiscales"]["layout"][2], "transform")
    @test schema_errors(s["multiscales"], bad) !== nothing
    bad = deepcopy(node)
    bad["attributes"]["multiscales"]["layout"][1]["asset"] = ""
    @test schema_errors(s["multiscales"], bad) !== nothing
end

@testitem "Matches the geospatial pyramid example of the multiscales convention" setup = [ConventionHelpers] begin
    using PyramidScheme: PyramidScheme as PS
    using DimensionalData
    using JSON
    using Statistics

    example = JSON.parsefile(joinpath(@__DIR__, "schemas", "examples", "geospatial-pyramid.json"))["attributes"]
    # 10 m cells covering the example bbox, north up
    dims = (X(500005.0:10.0:599995.0), Y(4999995.0:-10.0:4900005.0))
    factors = [1, 2, 4]
    attrs = group_node(PS.multiscales_group_attrs(pyramid_dims(dims, factors), factors;
        crs="EPSG:32632", resampling_method=mean))["attributes"]

    @test attrs["zarr_conventions"] == example["zarr_conventions"]
    # Every key of the example is reproduced, extra keys are allowed by the conventions
    for key in ("proj:code", "spatial:dimensions", "spatial:transform", "spatial:bbox")
        @test attrs[key] == example[key]
    end
    @test attrs["multiscales"]["resampling_method"] == example["multiscales"]["resampling_method"]
    for (ours, theirs) in zip(attrs["multiscales"]["layout"], example["multiscales"]["layout"])
        for (key, value) in theirs
            @test ours[key] == value
        end
    end
end

@testitem "Agrees with topozarr" setup = [ConventionHelpers] begin
    using PyramidScheme: PyramidScheme as PS
    using DimensionalData
    using JSON
    using Statistics

    fixtures = joinpath(@__DIR__, "fixtures", "topozarr")
    spatial_dims = (Dim{:x}, Dim{:y})
    # Same datasets as fixtures/topozarr/generate.py: north up, x and y starting at 0
    function griddims(nx, ny; time=nothing, res=1.0)
        dims = (Dim{:x}(range(0.0; step=res, length=nx)), Dim{:y}(range((ny - 1) * res; step=-res, length=ny)))
        time === nothing ? dims : (dims..., Dim{:time}(0:time-1))
    end
    cases = [
        ("square", griddims(32, 32), [1, 2, 4], "EPSG:4326", mean),
        ("non_square", griddims(48, 32), [1, 2, 4], "EPSG:4326", mean),
        ("extra_dim", griddims(32, 32; time=3), [1, 2, 4], "EPSG:4326", mean),
        ("projected", griddims(32, 32; time=2, res=100.0), [1, 2], "EPSG:3857", mean),
        ("sparse_factors", griddims(64, 48), [1, 2, 8], "EPSG:4326", mean ∘ skipmissing),
        ("sum", griddims(32, 32), [1, 2], "EPSG:4326", sum),
        ("odd_size", griddims(50, 30), [1, 2, 4], "EPSG:4326", mean),
    ]
    s = schemas()
    for (name, dims, factors, crs, method) in cases
        expected = JSON.parsefile(joinpath(fixtures, "$name.json"))
        # topozarr drops partial windows at the edges
        leveldims = [PS.level_dims(dims, f; edges=:drop, spatial_dims) for f in factors]
        attrs = PS.multiscales_group_attrs(leveldims, factors;
            crs, resampling_method=method, spatial_dims)
        node = group_node(attrs)
        ours = node["attributes"]
        @test sort(collect(keys(ours))) == sort(collect(keys(expected)))
        for key in setdiff(keys(expected), ["multiscales"])
            @test ours[key] == expected[key]
        end
        # topozarr names mean "mean" instead of "average" and leaves out the transform of level 0
        method_name(theirs) = get(Dict("mean" => "average"), theirs, theirs)
        @test ours["multiscales"]["resampling_method"] == method_name(expected["multiscales"]["resampling_method"])
        @test length(ours["multiscales"]["layout"]) == length(expected["multiscales"]["layout"])
        for (i, (level, theirs)) in enumerate(zip(ours["multiscales"]["layout"], expected["multiscales"]["layout"]))
            level, theirs = Dict{String,Any}(level), Dict{String,Any}(theirs)
            if i == 1
                transform = pop!(level, "transform")
                @test transform["scale"] == ones(length(dims))
                @test transform["translation"] == zeros(length(dims))
            else
                @test level["resampling_method"] == method_name(theirs["resampling_method"])
                level["resampling_method"] = theirs["resampling_method"]
            end
            @test level == theirs
        end
        # topozarr output is valid under the convention schemas as well
        for conv in ("multiscales", "proj", "spatial")
            @test schema_errors(s[conv], node) === nothing
        end
    end
end
