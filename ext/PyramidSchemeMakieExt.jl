module PyramidSchemeMakieExt
using Makie: Axis, Colorbar, DataAspect, Figure, FigureAxisPlot, Observable, Relative
using Makie: on, heatmap!, image!
import Makie: plot, plot!
using PyramidScheme: Pyramid, switchkeys, levels, selectlevel, xkey, ykey
using DimensionalData: DimensionalData as DD
using DimensionalData.Dimensions: XDim, YDim
using Extents: Extent, extent, intersects
miss2nan(x) = ismissing(x) ? NaN : x
"""
    plot(pyramids)
Plot a Pyramid. 
This will plot the coarsest resolution level at the beginning and will plot higher resolutions after zooming into the plot.
This is expected to be used with interactive Makie backends.
"""
function plot(pyramid::Pyramid;colorbar=true, size=(1155, 1043), kwargs...)
    #This should be converted into a proper recipe for Makie but this would depend on a pyramid type.
    fig = Figure(;size)
    lon, lat = DD.dims(DD.parent(pyramid))
    ax = Axis(fig[1,1], limits=(extrema(lon), extrema(lat)), aspect=DataAspect())
    hmap = plot!(ax, pyramid;kwargs...)
    if colorbar
        Colorbar(fig[1,2], hmap, height = Relative(3.5 / 4))
    end
    ax.autolimitaspect = 1
    FigureAxisPlot(fig, ax, hmap)
end

function plot!(ax, pyramid::Pyramid;interp=false, kwargs...)#; rastercrs=crs(parent(pyramid)),plotcrs=EPSG(3857), kwargs...)
    tip = levels(pyramid)[end-2][:,:]
    #@show typeof(tip)
    subtypes = Union{typeof.(pyramid.levels)..., typeof(parent(pyramid))}
    data = Observable{DD.AbstractDimMatrix}(tip)
    xval = only(values(extent(pyramid, XDim)))
    yval = only(values(extent(pyramid, YDim)))
    rasdataext = Extent(X=xval, Y=yval)
    rasext = extent(pyramid)
    on(ax.scene.viewport) do viewport
        limext = extent(ax.finallimits[])

        datalimit = switchkeys(limext, rasext)
        data.val = miss2nan.(selectlevel(pyramid, datalimit, target_imsize=viewport.widths))
        notify(data)
    end
    on(ax.finallimits) do limits
        limext = extent(limits)
        # Compute limit in raster projection
        #trans = Proj.Transformation(plotcrs, rastercrs, always_xy=true)
        #datalimit = trans_bounds(trans, limext)
        datalimit = switchkeys(limext, rasext)
        if intersects(rasdataext, limext)
            rasdata = selectlevel(pyramid, datalimit, target_imsize=ax.scene.viewport[].widths)
            # Project selected data to plotcrs
            #data.val = Rasters.resample(rasdata, crs=plotcrs, method=:bilinear )
            data.val = miss2nan.(rasdata)
        end
        notify(data)
    end
    #@show typeof(data)
    hmap = heatmap!(ax, data; interpolate=interp, kwargs...)
end
end