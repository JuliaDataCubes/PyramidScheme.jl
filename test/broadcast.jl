@testitem "broadcasting" begin
    using PyramidScheme: PyramidScheme as PS
    using DimensionalData
    using Rasters
    data = zeros(2000,2000)
    dd = DimArray(data, (X(11:2010), Y(101:2100)))
    pyramid = PS.Pyramid(dd)
    p0 = pyramid .- pyramid
    @test all(iszero, p0.base)
    @test all(all.(iszero, p0.levels))

    p1 = p0 .+ 1
    @test sum(p1.base) == prod(size(p1.base))
    for l in p1.levels
        @test sum(l) == prod(size(l))
    end
    tname = tempname() * ".tif"
    r = Raster(dd)
    write(tname, r, driver="cog", force=true)
    ptif = Pyramid(tname)
    # This fails because the pyramids have different layers
    @test_broken p0mix = ptif .- pyramid
 
end