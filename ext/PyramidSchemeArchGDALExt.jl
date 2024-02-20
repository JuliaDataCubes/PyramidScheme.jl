module PyramidSchemeArchGDALExt
using ArchGDAL: ArchGDAL as AG
using PyramidScheme: PyramidScheme as PS
using Rasters

function PS.Pyramid(path::AbstractString)
    base = Raster(path, lazy=true)
    agbase = AG.readraster(path)
    band = AG.getband(agbase,1)
    numlevels = AG.noverview(band)

    pyrlevels = if AG.nraster(agbase) > 1
        banddim = dims(base, 3)
        bandindices = 1: AG.nraster(agbase)
        dbase = dims(base,(1,2))
        pyrlevels = Raster[]
        for n in 0:numlevels-1
            levelbands = [Raster(AG.getoverview(AG.getband(agbase,bind), n), PS.agg_axis.(dbase, 2^(n+1))) for bind in bandindices]
            level = cat(levelbands..., dims=banddim)
            push!(pyrlevels, level)
        end
        pyrlevels
    else 
        pyrlevels = [Raster(AG.getoverview(band, n), PS.agg_axis.(dims(base), 2^(n+1))) for n in 0:numlevels-1]
    end
    PS.Pyramid(base, pyrlevels)
end

end