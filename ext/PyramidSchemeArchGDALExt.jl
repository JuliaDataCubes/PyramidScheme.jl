module PyramidSchemeArchGDALExt
using ArchGDAL: ArchGDAL as AG
using PyramidScheme: PyramidScheme as PS
import PyramidScheme: _pyramid_gdal
using YAXArrays
using DimensionalData


function _pyramid_gdal(path::AbstractString)
    base = Cube(path)
    agbase = AG.readraster(path)
    band = AG.getband(agbase,1)
    numlevels = AG.noverview(band)

    pyrlevels = if AG.nraster(agbase) > 1
        banddim = dims(base, 3)
        bandindices = 1: AG.nraster(agbase)
        dbase = dims(base,(1,2))
        pyrlevels = AbstractDimArray[]
        for n in 0:numlevels-1
            levelbands = [YAXArray(PS.agg_axis.(dbase, 2^(n+1)), AG.getoverview(AG.getband(agbase,bind), n)) for bind in bandindices]
            level = cat(levelbands..., dims=banddim)
            push!(pyrlevels, level)
        end
        pyrlevels
    else 
        pyrlevels = [YAXArray(PS.agg_axis.(dims(base), 2^(n+1)),AG.getoverview(band, n)) for n in 0:numlevels-1]
    end
    PS.Pyramid(base, pyrlevels, metadata(base))
end

end