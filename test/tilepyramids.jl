import PyramidScheme as PS
using Tyler: Map, wgs84
using Extents: Extent
import GLMakie
using HTTP
@dim lon XDim "Longitude"

@dim lat YDim "Latitude"
p2020 = PS.Pyramid("https://s3.bgc-jena.mpg.de:9000/pyramids/ESACCI-BIOMASS-L4-AGB-MERGED-100m-2020-fv4.0.zarr")
p2017 = PS.Pyramid("https://s3.bgc-jena.mpg.de:9000/pyramids/ESACCI-BIOMASS-L4-AGB-MERGED-100m-2017-fv4.0.zarr")
agbdiff = p2017 .- p2020
datarange = (0.0,400.0)

pv = PS.PyramidProvider(p2020,datarange...,colorscheme=:speed)

ext = Extent(X=(-180.0,180.0),Y=(-60.0,80.0))


s = HTTP.serve!(pv,"127.0.0.1",8765)


datarange = (-20. , 20.)
diffprov = PS.PyramidProvider(agbdiff, datarange..., colorscheme=:diverging_gkr_60_10_c40_n256)
s = HTTP.serve!(diffprov,"127.0.0.1",8766)
prov = TileProviders.Provider("http://127.0.0.1:8765/{z}/{x}/{y}.png")