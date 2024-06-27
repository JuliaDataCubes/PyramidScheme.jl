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
using PyramidScheme: PyramidScheme as PS
path = 
ras = Raster(EarthEnv{LandCover}, 6, lazy=true)
pyr = PS.Pyramid(ras)

# Now we can plot the data in GLMakie and get a nice interactive plot which uses the pyramids to provide a nice smooth experience by only loading the pixels which can fit into the Makie axis from an appropriate pyramid.

using GLMakie
plot(pyr)


```