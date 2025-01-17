@testitem "broadcasting" begin
    using PyramidScheme: PyramidScheme as PS
    using DimensionalData
    using Rasters
    #using ArchGDAL
    data = rand(2000,2000)
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
    #tname = tempname() * ".tif"
    #r = Raster(dd)
    #write(tname, r, driver="cog", force=true)
    #ptif = Pyramid(tname)
    #ptif0 = ptif .- ptif
    #@test all(iszero, ptif0.base)
    #@test all(all.(iszero, ptif0.levels))

    #pyr500 = Pyramid(dd)
    # This fails because the pyramids have different layers
    #p0mix = ptif .- pyr500
    #@test all(iszero, p0mix.base)
 
end