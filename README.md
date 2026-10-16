# PyramidScheme
[![codecov](https://codecov.io/gh/JuliaDataCubes/PyramidScheme.jl/graph/badge.svg?token=EBhwEBXlnx)](https://codecov.io/gh/JuliaDataCubes/PyramidScheme.jl)


PyramidScheme.jl is a package to easily and efficiently compute pyramids of a given datacube which might be larger than RAM.
It uses DiskArrayEngine.jl as the computational backend.
The pyramids can then be used to interactively explore the data.

The long term aim of PyramidScheme.jl is to enable computing based on the layers of the pyramid to enable a more interactive exploration of computations based of pyramided datasets. 

[PyramidScheme.jl in action](https://github.com/JuliaDataCubes/PyramidScheme.jl/assets/17124431/63166348-7b18-4af1-b6bb-4eb0af2def9e)


## Usage

To compute the pyramids of a given data cube use the following steps:

```julia
using Rasters
using RasterDataSources
using ArchGDAL
using Statistics
using Extents
using PyramidScheme: PyramidScheme as PS
ras = Raster(WorldClim{Elevation},:elev, res="30s", lazy=true)
pyr = PS.Pyramid(ras)

# Now we can plot the data in GLMakie and get a nice interactive plot which uses the pyramids to provide a nice smooth experience by only loading the pixels which can fit into the Makie axis from an appropriate pyramid.

using GLMakie
plot(pyr)
```

## Writing pyramids to Zarr

`buildpyramids(src, dest)` computes the pyramid levels of `src` and writes them as a multiscale Zarr v3 store at `dest`.
`src` can be the path of a dataset that YAXArrays can open or any `AbstractDimArray`, and it is not modified.
The store follows the [multiscales](https://github.com/zarr-conventions/multiscales), [spatial](https://github.com/zarr-conventions/spatial) and [proj](https://github.com/zarr-conventions/proj) conventions:
the root group only holds the metadata, `src` is copied to the child group `0` and the aggregated levels are written to the groups `1`, `2`, ...

Here we use the same example data as above:

```julia
using YAXArrays
elevpath = getraster(WorldClim{Elevation},:elev, res="30s")
elev = Cube(elevpath)
pyrpath = tempname() * ".zarr"
@time pyr = PS.buildpyramids(elev, pyrpath; resampling_method=mean, crs="EPSG:4326")
```

The most important keywords are:

- `nlevels` or `factors`: the levels to build, e.g. `nlevels=3` for factors `[1, 2, 4, 8]` or `factors=[1, 4, 16]`.
  By default enough levels are built for the coarsest one to fit into 256 by 256 cells.
- `resampling_method`: the function aggregating the cells of a window, e.g. `mean`, `sum` or `mean ∘ skipmissing`.
- `edges`: `:keep` (default) turns a partial window at the end of a dimension into its own cell, `:drop` leaves it out.
- `crs`: the coordinate reference system written to the `proj:` attributes. By default it is taken from the dimensions, e.g. of a Rasters.jl `Raster`.

An in memory `Pyramid` can be written the same way with `write(dest, pyramid)`.

`Pyramid(path)` opens a multiscale store, written by PyramidScheme or by other tools such as topozarr, as well as stores written by earlier versions of PyramidScheme.
Remote stores need the corresponding storage package, e.g. `using Zarr` loads all of them.

### Example above

To reproduce the zooming example shown at the top of the README you can get the Above Ground Biomass here or try to load it from the cloud, but this might have some serious lag.
```julia
using YAXArrays, PyramidScheme
using Zarr # for reading from http
using GLMakie
p2020 = PS.Pyramid("https://s3.bgc-jena.mpg.de:9000/pyramids/ESACCI-BIOMASS-L4-AGB-MERGED-100m-2020-fv4.0.zarr")
replacenan(nanval) =  data -> <=(nanval)(data) ? NaN32 : Float32(data)
p2020nan = replacenan(0).(p2020)
plot(p2020nan, colormap=:speed, colorscale=sqrt)

# to save the base data use YAXArrays to load the underlying data
c = Cube("https://s3.bgc-jena.mpg.de:9000/pyramids/ESACCI-BIOMASS-L4-AGB-MERGED-100m-2020-fv4.0.zarr")
savecube(c, "somepath.zarr")
# Build the pyramid from this data locally
buildpyramids("somepath.zarr", "pyramid.zarr")
```
