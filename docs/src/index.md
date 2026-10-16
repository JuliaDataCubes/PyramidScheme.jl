# PyramidScheme.jl

PyramidScheme.jl computes pyramids of data cubes that might be larger than memory and uses them to
explore the data interactively. A [`Pyramid`](@ref) behaves like the array of its highest resolution,
subsetting it returns another `Pyramid`, and plotting it with Makie only loads the level that fits
into the axis.

```julia
using PyramidScheme, DimensionalData
pyr = Pyramid(DimArray(rand(4000, 3000), (X(1:4000), Y(1:3000))))
```

Pyramids can be written to and read from Zarr stores that follow the multiscales convention,
see [Multiscale Zarr stores](@ref), and Cloud Optimized GeoTIFFs with overviews can be opened
with `Pyramid(path)` after loading ArchGDAL.

```@docs
Pyramid
```
