@testitem "Level grids" begin
    using PyramidScheme: PyramidScheme as PS
    using DimensionalData
    using DimensionalData.Lookups

    x = X(0.0:9.0)
    @test collect(PS.coarsen_dim(x, 4)) == [1.5, 5.5, 9.5]
    @test collect(PS.coarsen_dim(x, 4; edges=:drop)) == [1.5, 5.5]
    @test PS.coarsen_dim(x, 1) === x
    y = Y(9.0:-1.0:0.0)
    @test collect(PS.coarsen_dim(y, 4)) == [7.5, 3.5, -0.5]
    @test collect(PS.coarsen_dim(y, 4; edges=:drop)) == [7.5, 3.5]

    # Levels keep the origin of the base and have factor times larger cells
    base = dims(DimArray(rand(10, 7), (X(0.0:9.0), Y(6.0:-1.0:0.0))))
    t0 = PS.spatial_transform(base...)
    @test_throws ArgumentError PS.level_dims(base, 8; edges=:drop)
    for f in (2, 4, 8), edges in (:keep, :drop)
        f == 8 && edges === :drop && continue
        lx, ly = PS.level_dims(base, f; edges)
        @test length(lx) == (edges === :keep ? cld(10, f) : 10 ÷ f)
        @test length(ly) == (edges === :keep ? cld(7, f) : 7 ÷ f)
        t = PS.spatial_transform(lx, ly; fallback_res=(t0[1] * f, t0[5] * f))
        @test t ≈ [t0[1] * f, 0.0, t0[3], 0.0, t0[5] * f, t0[6]]
        # The lookup keeps its type and order
        @test lookup(ly) isa Sampled
        @test order(ly) isa ReverseOrdered
        @test span(lx) == Regular(f * 1.0)
    end

    # Interval lookups keep their locus
    start = X(Sampled(0.0:10.0:90.0; sampling=Intervals(Start()), order=ForwardOrdered(), span=Regular(10.0)))
    cstart = PS.coarsen_dim(start, 2)
    @test collect(cstart) == 0.0:20.0:80.0
    @test locus(cstart) isa Start
    @test span(cstart) == Regular(20.0)
    stop = X(Sampled(10.0:10.0:100.0; sampling=Intervals(End()), order=ForwardOrdered(), span=Regular(10.0)))
    @test collect(PS.coarsen_dim(stop, 2)) == 20.0:20.0:100.0

    @test_throws ArgumentError PS.coarsen_dim(X(0.0:2.0), 4; edges=:drop)
    @test_throws ArgumentError PS.coarsen_dim(x, 2; edges=:trim)

    # Only spatial dimensions are downsampled
    withtime = (X(0.0:9.0), Y(0.0:9.0), Ti(1:3))
    @test length.(PS.level_dims(withtime, 4)) == (3, 3, 3)
    @test PS.level_dims(withtime, 4)[3] == withtime[3]
    @test length.(PS.level_dims((Dim{:lon}(0.0:9.0), Dim{:lat}(0.0:9.0)), 2;
        spatial_dims=(Dim{:lon}, Dim{:lat}))) == (5, 5)
end

@testitem "Pyramid factors and level indices" begin
    using PyramidScheme: PyramidScheme as PS
    using DimensionalData

    @test PS.pyramid_factors() == [1]
    @test PS.pyramid_factors(default_nlevels=3) == [1, 2, 4, 8]
    @test PS.pyramid_factors(nlevels=2, default_nlevels=5) == [1, 2, 4]
    @test PS.pyramid_factors(factors=[1, 2, 8]) == [1, 2, 8]
    @test_throws ArgumentError PS.pyramid_factors(nlevels=2, factors=[1, 2])
    @test_throws ArgumentError PS.pyramid_factors(nlevels=-1)
    @test_throws ArgumentError PS.pyramid_factors(factors=[1, 3, 4])

    @test PS.levelindex(4, 9, 3) == 3
    @test PS.levelindex(4, 12, 2) == 2
    @test PS.levelindex(4, 1:8, 3) == 1:2
    @test PS.levelindex(4, :, 3) == Colon()

    x = X(0.0:9.0)
    @test PS.level_factor(x, PS.coarsen_dim(x, 4)) == 4
    @test PS.level_factor(x, PS.coarsen_dim(x, 4; edges=:drop)) == 4
    @test PS.level_factor(x, x) == 1
    @test PS.level_factor(x, PS.coarsen_dim(x, 16)) == 10
end

@testitem "Aggregation with partial windows" begin
    using PyramidScheme: PyramidScheme as PS
    using DimensionalData
    using Statistics

    ones2d = DimArray(ones(10, 6), (X(0.0:9.0), Y(0.0:5.0)))
    keep = Pyramid(ones2d; resampling_method=sum, factors=[1, 4])
    @test PS.nlevels(keep) == 1
    @test parent(keep.levels[1]) == [16 8; 16 8; 8 4]
    drop = Pyramid(ones2d; resampling_method=sum, factors=[1, 4], edges=:drop)
    @test parent(drop.levels[1]) == fill(16.0, 2, 1)
    @test collect(dims(drop.levels[1], X)) == [1.5, 5.5]

    # Recursive and direct aggregation agree for sum
    rec = Pyramid(ones2d; resampling_method=sum, factors=[1, 2, 4])
    direct = Pyramid(ones2d; resampling_method=sum, factors=[1, 2, 4], recursive=false)
    @test parent(rec.levels[2]) == parent(direct.levels[2]) == [16 8; 16 8; 8 4]

    data = rand(10, 6)
    dd = DimArray(data, (X(0.0:9.0), Y(0.0:5.0)))
    for recursive in (true, false)
        pyr = Pyramid(dd; resampling_method=mean, factors=[1, 2, 8], recursive)
        level = parent(pyr.levels[2])
        @test size(level) == (2, 1)
        @test level[1, 1] ≈ mean(data[1:8, 1:6])
        @test level[2, 1] ≈ mean(data[9:10, 1:6])
    end
    # Y has 6 cells, too few for a factor 8 level without partial windows
    @test_throws ArgumentError Pyramid(dd; resampling_method=mean, factors=[1, 2, 8], edges=:drop)
    dropped = Pyramid(dd; resampling_method=mean, factors=[1, 2, 4], edges=:drop)
    @test size.(dropped.levels) == [(5, 3), (2, 1)]
    @test parent(dropped.levels[2])[2, 1] ≈ mean(data[5:8, 1:4])
end

@testitem "Pyramids with an extra dimension" begin
    using PyramidScheme: PyramidScheme as PS
    using DimensionalData

    data = ones(16, 8, 3)
    for t in 1:3
        data[:, :, t] .= t
    end
    dd = DimArray(data, (X(0.0:15.0), Y(0.0:7.0), Ti(1:3)))
    pyr = Pyramid(dd; resampling_method=sum, factors=[1, 2, 4])
    @test size.(pyr.levels) == [(8, 4, 3), (4, 2, 3)]
    @test all(parent(pyr.levels[1])[:, :, 2] .== 8)
    @test all(parent(pyr.levels[2])[:, :, 3] .== 48)

    # Selecting along the extra dimension selects the same step in every level
    sub = pyr[Ti=2]
    @test sub isa Pyramid
    @test size.(sub.levels) == [(8, 4), (4, 2)]
    @test all(parent(sub.levels[1]) .== 8)
    @test all(parent(sub.levels[2]) .== 32)

    window = pyr[1:8, 1:4, 3]
    @test size.(window.levels) == [(4, 2), (2, 1)]
    @test all(parent(window.levels[2]) .== 48)
end

@testitem "Subsetting pyramids with sparse factors and dropped edges" begin
    using PyramidScheme: PyramidScheme as PS
    using DimensionalData

    pyr = Pyramid(DimArray(rand(64, 48), (X(0.0:63.0), Y(0.0:47.0))); factors=[1, 2, 8])
    sub = pyr[1:16, 1:16]
    @test size.(sub.levels) == [(8, 8), (2, 2)]

    dropped = Pyramid(DimArray(rand(10, 10), (X(0.0:9.0), Y(0.0:9.0))); factors=[1, 4], edges=:drop)
    @test size(dropped.levels[1]) == (2, 2)
    # Indices in the dropped edge map to the last cell of the level
    edge = dropped[9:10, 1:4]
    @test size(edge.levels[1]) == (1, 1)
end
