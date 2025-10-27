import PyramidScheme as PS
#using Tyler: Map, wgs84
using Extents: Extent
#import GLMakie
using HTTP
using DimensionalData.Dimensions
@dim lon XDim "Longitude"

@dim lat YDim "Latitude"
replacenan(data) =  data <= 0 ? NaN32 : Float32(data)
p2020 = replacenan.(PS.Pyramid("https://s3.bgc-jena.mpg.de:9000/pyramids/ESACCI-BIOMASS-L4-AGB-MERGED-100m-2020-fv4.0.zarr"))
p2017 = replacenan.(PS.Pyramid("https://s3.bgc-jena.mpg.de:9000/pyramids/ESACCI-BIOMASS-L4-AGB-MERGED-100m-2017-fv4.0.zarr"))
agbdiff = p2017 .- p2020
datarange = (1.0,200.0)

pv = PS.PyramidProvider(PS.cache(p2020),datarange...,colorscheme=:speed)

ext = Extent(X=(-180.0,180.0),Y=(-60.0,80.0))


s = HTTP.serve!(pv,"127.0.0.1",8766)


datarange = (-15. , 15.)
diffprov = PS.PyramidProvider(PS.cache(agbdiff), datarange..., colorscheme=:diverging_gkr_60_10_c40_n256)
portnumber = rand(8700:8800)
sd = HTTP.serve!(diffprov,"127.0.0.1",portnumber)
prov = TileProviders.Provider("http://127.0.0.1:$portnumber/{z}/{x}/{y}.png")