using Test, TestItemRunner

@run_package_tests@testset "Pyramid" begin
    using DimensionalData
    using PyramidScheme: PyramidScheme as PS
    using Makie
    data = rand(2000,2000)
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
    @test PS.compute_nlevels(zeros(1000)) == 0
    @test PS.compute_nlevels(zeros(1000,1025)) == 1
    @test PS.compute_nlevels(zeros(10000, 8000)) == 4
end

@testitem "Broadcast" begin
    using DimensionalData
    using PyramidScheme: PyramidScheme as PS
    data1 = rand(2000,2000)
    dd1 = DimArray(data1, (X(1:2000), Y(1:2000)))
    pyramid1 = PS.Pyramid(dd1)
    data2 = rand(2000,2000)
    dd2 = DimArray(data2, (X(1:2000), Y(1:2000)))
    pyramid2 = PS.Pyramid(dd2)
    pyrsum = pyramid1 .+ pyramid2
    @test pyrsum.levels[end][1,1] == pyramid1.levels[end][1,1] + pyramid2.levels[end][1,1]
end
