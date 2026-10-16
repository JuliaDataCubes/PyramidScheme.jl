
File filter
# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Links are updated with Changelog.jl using the command:

```julia
Changelog.generate(
    Changelog.CommonMark(),             # output type
    "CHANGELOG.md";                     # input and output file
    repo = "JuliaDataCubes/PyramidScheme.jl", # default repository for links
)
```

## [0.2.0]

### Breaking

- `buildpyramids(src, dest)` writes the pyramid to a new store at `dest`. Building the levels inside the
  source store, `buildpyramids(path)`, throws an error, because the multiscales convention stores every
  level, including the base, in its own child group.
- Pyramids are written as Zarr v3 stores following the multiscales, spatial and proj conventions
  (`v0.1` of https://github.com/zarr-conventions) instead of the previous `multiscales` attribute.
- Depend on ZarrCore, ZarrZstd and ZarrBlosc instead of Zarr. Reading from remote stores needs the
  corresponding storage package, e.g. `using Zarr`.
- The metadata of `Pyramid(data)` records the convention name of the resampling method, e.g.
  `"average"` for `mean`.
- `buildpyramids` aggregates with `mean ∘ skipmissing` by default instead of `mean`, like `Pyramid(data)`
  and topozarr, so missing values are left out instead of making the coarse cell missing.

### Added

- `nlevels` and `factors` keywords to choose the levels, e.g. `factors=[1, 4, 16]`.
- `edges=:drop` to leave out partial windows at the end of a dimension instead of keeping them as a cell.
- `crs` keyword for the `proj:` attributes, taken from the dimensions by default.
  They are written to the root group and to every level group.
- `write(dest, pyramid)` to write an in memory pyramid.
- `Pyramid(path)` reads stores with a multiscales layout, including stores written by topozarr.
  Stores written by earlier versions can still be read.
- Consolidated metadata and zstd compression for written stores.

### Fixed

- The coordinates of the levels lie exactly on the grid of the level, also when the length of a
  dimension is not a multiple of the factor.
- Subsetting a pyramid along a non spatial dimension, e.g. `pyr[Ti=2]`, selects the same step in all
  levels.
- In memory pyramids only aggregate the spatial dimensions.
- `selectlevel` works with factors that are not powers of two.

## [0.1.2]

### Added 

- Implemented proper MakieRecipies this should make Observable(::Pyramid) work
- Added `cat` of Pyramids to concatenate on all levels

### Fixed

- Improved computation of the dimensions of the aggregated levels
