@testitem "broadcasting two Arrays" begin
    using PyramidScheme: PyramidScheme as PS
    using DimensionalData
    using ArchGDAL

    pyr = Pyramid("test/data/pyramidmiddle.tif")
    p0 = pyr .- pyr
    @test all(iszero, p0.base)
    @test all(all.(iszero, p0.levels))

    p1 = p0 .+ 1
    @test sum(p1.base) == prod(size(p1.base))
    for l in p1.levels
        @test sum(l) == prod(size(l))
    end
    
end