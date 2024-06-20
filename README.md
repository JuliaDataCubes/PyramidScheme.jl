# PyramidScheme

PyramidScheme.jl is a package to easily and efficiently compute pyramids of a given datacube which might be larger than RAM.
It uses DiskArrayEngine.jl as the computational backend.
The pyramids can then be used to interactively explore the data.

The long term aim of PyramidScheme.jl is to enable computing based on the layers of the pyramid to enable a more interactive exploration of computations based of pyramided datasets. 



## Usage

To compute the pyramids of a given data cube use the following steps:

```julia
using Rasters
using RasterDataSources
using ArchGDAL
using Statistics
using Extents
using Tyler
using PyramidScheme: PyramidScheme as PS
ras = Raster(WorldClim{Weather}, :tmax, date=Date(2018, 6,1))

n_level = PS.compute_nlevels(ras)

# Compute the size of every pyramid and the axes of the pyramids.
pyramid_sizes =  [ceil.(Int, size(ras) ./ 2^i) for i in 1:n_level]
pyramid_axes = [PS.agg_axis.(dims(ras),2^i) for i in 1:n_level]

# Define the output arrays, these will be saved to disk in zarr format if they are too large. 
outarrs = PS.output_arrays(pyramid_sizes, Float32)

# Now we can compute the pyramids and put them into outarrs

@time PS.fill_pyramids(ras.data,outarrs,mean,true)

# And we put the pyramids and the original data into one vector for later access

pyramids = [ras, Raster.(outarrs, pyramid_axes)...]


# Now we can plot the data in GLMakie and get a nice interactive plot which uses the pyramids to provide a nice smooth experience by only loading at most 1024 by 512 pixels from an appropriate pyramid.

using GLMakie

lon,lat = dims(ras)

fig = Figure()

ax = Axis(fig[1,1], limits=(extrema(lon), extrema(lat)))
data = Observable{Raster}(pyramids[end])

on(ax.finallimits) do limits
    data.val = PS.selectlevel(pyramids, extent(limits), lon.val.span.step)
    @show size(data.val)
    notify(data)
end

heatmap!(ax, data, interpolate=false)
ax.autolimitaspect=1

```