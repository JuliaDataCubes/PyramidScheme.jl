using PyramidScheme
using GLMakie
using DimensionalData.Dimensions

@dim lat YDim "Latitude"
@dim lon XDim "Longitude"
p1 = Pyramid("test/data/ESACCI-BIOMASS-L4-AGB-MERGED-100m-2020-fv4.0.zarr/")
p2 = Pyramid("test/data/ESACCI-BIOMASS-L4-AGB-MERGED-100m-2017-fv4.0.zarr/")
agbdiff = p1 .- p2
#agbdiff[agbdiff.==0] .= NaN

f, a, p = plot(agbdiff, colormap=Reverse(:diverging_gwr_55_95_c38_n256))
lift(p.converted[3]) do ras
    actual_crange = Makie.distinct_extrema_nan(ras)
    maxval = maximum(abs.(actual_crange))
    p.colorrange[] = (-maxval, maxval)
end

agbrel = (abs.(agbdiff) .> 10) .* agbdiff
f, a, p = plot(agbrel, colormap=Reverse(:diverging_gwr_55_95_c38_n256))
lift(p.converted[3]) do ras
    actual_crange = Makie.distinct_extrema_nan(ras)
    maxval = maximum(abs.(actual_crange))
    p.colorrange[] = (-maxval, maxval)
end