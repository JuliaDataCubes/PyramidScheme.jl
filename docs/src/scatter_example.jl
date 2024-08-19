using GLMakie
using DimensionalData.Dimensions
using PyramidScheme
using Extents
using Rasters
using RasterDataSources
using ArchGDAL
@dim lat YDim "Latitude"
@dim lon XDim "Longitude"
replacenan(data) =  data <= 0 ? NaN32 : Float32(data)


rastas = Raster(CHELSA{Climate},:tas, lazy=true, month=1)
rastasscale = rastas.metadata["scale"] .* rastas .+ rastas.metadata["offset"]
pyrtas = Pyramid(rastasscale)

rasprecip = Raster(CHELSA{Climate},:pr, lazy=true, month=1)
rasprecipscale = rasprecip.metadata["scale"] .* rasprecip .+ rasprecip.metadata["offset"]

pyrprecip = Pyramid(rasprecip)

scattas = Observable(vec(collect(pyrtas.levels[end])))
scatprecip = Observable(vec(collect(pyrprecip.levels[end])))

scatpoints = Observable(Point2f[])

begin
figmap = Figure()
axtas = Axis(figmap[1,1])
plot!(axtas, pyrtas, colormap=:berlin);
axprecip = Axis(figmap[1,2])
plot!(axprecip, pyrprecip)
axscat = Axis(figmap[2,1:2])
#scatter!(axscat, scatpoints, scatprecip)
#datashader!(axscat, scatpoints)
on(axtas.finallimits) do limits
    limext = Extents.extent(limits)
    rasdataext = Extents.extent(pyrtas)
    datalimit = PyramidScheme.switchkeys(limext, rasdataext)
    if Extents.intersects(rasdataext, datalimit)
        scatpoints.val = Point2f.(vec(collect(PyramidScheme.selectlevel(pyrtas, datalimit, target_imsize=axtas.scene.viewport[].widths))),
                            vec(collect(PyramidScheme.selectlevel(pyrprecip, datalimit, target_imsize=axtas.scene.viewport[].widths))))

    end
    notify(scatpoints)
#    notify(scatprecip)
end
datashader!(axscat, scatpoints)
linkaxes!(axtas, axprecip)
figmap
end