module PyramidSchemeMakieExt
using Makie: Axis, Colorbar, DataAspect, Figure, FigureAxisPlot, Observable, Relative
using Makie: on, heatmap!, map!
using Makie.ComputePipeline: add_input!
using Makie

using PyramidScheme: Pyramid, switchkeys, levels, selectlevel, xkey, ykey
using DimensionalData: DimensionalData as DD
using DimensionalData.Dimensions: XDim, YDim
using Extents: Extent, extent, intersects
miss2nan(x) = ismissing(x) ? NaN : x

# hacks to get around DD hacks that get around Makie issues
for p in (Heatmap, Image, Contour, Contourf, Contour3d, Spy, Surface) 
    f = Makie.plotkey(p)
    @eval begin
        function Makie.$f(A::Pyramid; kwargs...)
            invoke(Makie.$f, Tuple{AbstractMatrix{<: Any}}, A; kwargs...)
        end
        function Makie.$f(A::Observable{<: Pyramid}; kwargs...)
            invoke(Makie.$f, Tuple{<:Union{<:AbstractMatrix{<: Any}, <:Observable{<: AbstractMatrix{<:Any}}}}, A; kwargs...)
        end
        Makie.expand_dimensions(::Type{$p}, p::Pyramid) = (p,)
        Makie.convert_arguments(::Type{$p}, p::Pyramid) = (p,)
    end
end

struct PyramidConversion <: Makie.ConversionTrait end

function Makie.conversion_trait(::Type{<:Heatmap}, ::Pyramid)
    return PyramidConversion()
end

function Makie.types_for_plot_arguments(::Type{<:Heatmap}, ::PyramidConversion)
    return Tuple{<:Pyramid}
end

Makie.expand_dimensions(::PyramidConversion, p::Pyramid) = (p,)
Makie.convert_arguments(::PyramidConversion, p::Pyramid) = (p,)



Makie.plottype(::Pyramid) = Makie.Heatmap

function Makie.plot!(plot::Heatmap{<: Tuple{<: Pyramid}})
    #=
    Go from a relative space rectangle 
    to the rectangle in data space and pixel space,
    thus getting `ax.finallimits` and `ax.viewport`
    respectively.

    This is a bit arcane but probably the simplest thing in the long run.
    =#
    inputpositions = [Point2f(0, 0), Point2f(1, 1)]
    add_input!(plot.attributes, :__pyramid_input_positions, inputpositions)
    Makie.register_positions_projected!(
        plot; input_space = :relative, output_space = :space,
        input_name = :__pyramid_input_positions,
        output_name = :__pyramid_dataspace_positions,
    )
    Makie.register_positions_projected!(
        plot; input_space = :relative, output_space = :pixel,
        input_name = :__pyramid_input_positions,
        output_name = :__pyramid_pixelspace_positions,
    )

    data = Observable{DD.AbstractDimMatrix}(levels(plot.arg1[])[end-2][:, :])
    onany(plot, plot.arg1, plot.__pyramid_dataspace_positions, plot.__pyramid_pixelspace_positions) do pyramid, datapos, pixelpos
        xyext = values.(extent(pyramid, (XDim, YDim)))
        xval, yval = first(xyext), last(xyext)
        pyramid_data_ext = Extent(X=xval, Y=yval)
        pyramid_ext = extent(pyramid)

        data_limits_ext = Extent(X = extrema(first, datapos), Y = extrema(x -> x[2], datapos))
        pixel_widths = Point2f(abs.(pixelpos[2] .- pixelpos[1]))

        datalimit = switchkeys(data_limits_ext, pyramid_ext)

        if intersects(pyramid_data_ext, data_limits_ext)
            data[] = miss2nan.(
                selectlevel(pyramid, datalimit, target_imsize = pixel_widths)
            )
        end
        nothing
    end

    heatmap!(plot, plot.attributes, data)
end

function Makie.data_limits(p::Heatmap{<: Tuple{<: Pyramid}})
    extX, extY = extent(p.arg1[], (XDim, YDim))
    rect = Rect3f((extX[1], extY[1],0), (extX[2] - extX[1], extY[2] - extY[1], 0))
    return rect
end
Makie.boundingbox(p::Heatmap{<: Tuple{<: Pyramid}}, space::Symbol = :data) = Makie.apply_transform_and_model(p, Makie.data_limits(p))


# """
#     plot(pyramids)
# Plot a Pyramid. 
# This will plot the coarsest resolution level at the beginning and will plot higher resolutions after zooming into the plot.
# This is expected to be used with interactive Makie backends.
# """
# function plot(pyramid::Pyramid;colorbar=true, size=(1155, 1043), kwargs...)
#     #This should be converted into a proper recipe for Makie but this would depend on a pyramid type.
#     fig = Figure(;size)
#     lon, lat = DD.dims(DD.parent(pyramid))
#     ax = Axis(fig[1,1], limits=(extrema(lon), extrema(lat)), aspect=DataAspect())
#     hmap = plot!(ax, pyramid;kwargs...)
#     if colorbar
#         Colorbar(fig[1,2], hmap, height = Relative(3.5 / 4))
#     end
#     ax.autolimitaspect = 1
#     FigureAxisPlot(fig, ax, hmap)
# end

# function plot!(ax, pyramid::Pyramid;interp=false, kwargs...)#; rastercrs=crs(parent(pyramid)),plotcrs=EPSG(3857), kwargs...)
#     tip = levels(pyramid)[end-2][:,:]
#     #@show typeof(tip)
#     subtypes = Union{typeof.(pyramid.levels)..., typeof(parent(pyramid))}
#     data = Observable{DD.AbstractDimMatrix}(tip)
#     xyext = values.(extent(pyramid, (XDim, YDim)))
#     xval = first(xyext)
#     yval = last(xyext)
#     rasdataext = Extent(X=xval, Y=yval)
#     rasext = extent(pyramid)
#     on(ax.scene.viewport) do viewport
#         limext = extent(ax.finallimits[])

#         datalimit = switchkeys(limext, rasext)
#         data[] = miss2nan.(selectlevel(pyramid, datalimit, target_imsize=viewport.widths))
#     end
#     on(ax.finallimits) do limits
#         limext = extent(limits)
#         # Compute limit in raster projection
#         #trans = Proj.Transformation(plotcrs, rastercrs, always_xy=true)
#         #datalimit = trans_bounds(trans, limext)
#         datalimit = switchkeys(limext, rasext)
#         if intersects(rasdataext, limext)
#             rasdata = selectlevel(pyramid, datalimit, target_imsize=ax.scene.viewport[].widths)
#             # Project selected data to plotcrs
#             #data.val = Rasters.resample(rasdata, crs=plotcrs, method=:bilinear )
#             data.val = miss2nan.(rasdata)
#         end
#         notify(data)
#     end
#     #@show typeof(data)
#     hmap = heatmap!(ax, data; interpolate=interp, kwargs...)
# end
end