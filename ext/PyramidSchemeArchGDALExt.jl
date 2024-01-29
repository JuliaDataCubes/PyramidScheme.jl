module PyramidSchemeArchGDALExt
using ArchGDAL: ArchGDAL as AG
using PyramdScheme: PyramidScheme as PS
using Rasters

function Pyramid(path::AbstractString)
    base = Raster(path, lazy=true)
    agbase = AG.readraster(path)
    band = AG.getband(agbase,1)
    numlevels = AG.noverview(band)
    dbase = dims(base)
    pyrlevels = [Raster(AG.getoverview(band, n), PS.agg_axis.(dims(base), 2^(n+1))) for n in 0:numlevels-1]
    PS.Pyramid(base, levels)
end



end