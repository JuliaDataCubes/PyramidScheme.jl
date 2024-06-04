# Aim of PyramidScheme.jl

The aim of the PyramidScheme.jl package is to provide an easy way to interactively work with larger than memory data.
When the data does not have pyramid layers computed it provides an easy way to compute the pyramid layers and save them efficiently in Zarr format. 
When the pyramid layers are already available this information will be used to run the computations first on the coarsest level to then enable a first glimpse at the results and to interactively run the computations on the parts of the data that the user is currently looking at. 

## Simple example

As a proof of concept we could look at the sum of two datasets which are too large to fit into memory and therefore the computation of the sum should be taking some time but when we only do it on the highest level which would by convention be only 256 by 256 pixels it should be able to be plotted instantaneously. 

```julia
using PyramidScheme
pyr = Pyramid("path/to/layeredfile.tiff")
pyrzarr = Pyramid("path/to/layeredfile.zarr")

# This computation should only be lazily registered at this stage
pyrsum = pyr .+ pyrzarr

using GLMakie
plot(pyrsum) # This should plot the sum in the highest level and zooming in would compute the rest of the layers.
```