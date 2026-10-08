# Multiscale Zarr stores

PyramidScheme writes pyramids as Zarr v3 stores that follow the
[multiscales](https://github.com/zarr-conventions/multiscales),
[spatial](https://github.com/zarr-conventions/spatial) and
[proj](https://github.com/zarr-conventions/proj) conventions, version `v0.1`.
These stores can be read by other tools that implement the conventions,
for example [topozarr](https://github.com/carbonplan/topozarr) and
[zarr-layer](https://github.com/carbonplan/zarr-layer).

## Building a pyramid

```julia
using PyramidScheme, YAXArrays, Statistics
pyr = buildpyramids("input.zarr", "pyramid.zarr"; resampling_method=mean, crs="EPSG:4326")
```

`buildpyramids(src, dest)` reads `src`, which is either the path of a dataset that YAXArrays can open
or an `AbstractDimArray`, and writes the store to `dest`. The source is not modified.
The result is the [`Pyramid`](@ref) opened from `dest`.

A pyramid that has been computed in memory, e.g. with `Pyramid(data)`, is written with
`write(dest, pyramid)`, which takes the same keywords.

## Layout of the store

```
pyramid.zarr/
├── zarr.json      # multiscales, spatial: and proj: attributes, consolidated metadata
├── 0/             # the source data
│   ├── layer
│   ├── X
│   └── Y
├── 1/             # downsampled by the second factor
│   └── ...
└── 2/
    └── ...
```

The root group only holds metadata. Every level, including the base, is stored in a child group named
after its index, as the multiscales convention requires.
The arrays carry `dimension_names` and the data arrays are compressed with zstd.

The root group has the following attributes:

- `zarr_conventions`: the conventions used by the store.
- `multiscales`: the `layout` of the levels. Every level after the base records the level it is
  `derived_from`, the `transform` (scale and translation) relative to that level,
  the `resampling_method`, and its `spatial:transform` and `spatial:shape`.
- `spatial:dimensions`, `spatial:transform`, `spatial:bbox`, `spatial:shape` and `spatial:registration`
  of the base.
- `proj:code` and `proj:wkt2`, when the CRS is known. `proj:code` is only written when the CRS has an
  authority code.

## Levels

The levels are chosen with `nlevels` or `factors`:

- `nlevels=3` builds three aggregated levels with the factors `[1, 2, 4, 8]`.
- `factors=[1, 4, 16]` builds levels downsampled by 4 and 16. Every factor must be a multiple of the
  previous one.
- Without either, enough levels are built for the coarsest one to fit into `tilesize` (256) cells.

With `recursive=true` (the default) every level is computed from the previous one, which is faster but
not exact for aggregations such as the median.

When the length of a spatial dimension is not a multiple of the factor, the last window is partial.
`edges=:keep` (the default) turns it into a cell of its own, so no data is lost at coarse levels.
`edges=:drop` leaves it out, so that every cell aggregates a full window.

Every level starts at the same corner as the base and has cells that are `factor` times larger,
so the coordinates of the levels lie exactly on their grid.

## Metadata styles

The conventions leave some room for interpretation, so the `multiscales` metadata can be written in
two styles:

| | `ConventionStyle()` (default) | `CellCenterStyle()` |
|---|---|---|
| `translation` between levels | `0`, as in the examples of the convention | `(step - 1) / 2`, measured between cell centers |
| `transform` of the base level | identity | left out |
| name of `mean` | `"average"` | `"mean"` |

`style=CellCenterStyle(), edges=:drop` writes the same metadata as topozarr.

## Coordinate reference system

By default (`crs=:auto`) the CRS is taken from the lookups of the spatial dimensions, for example of a
`Raster` from Rasters.jl. It can be given explicitly with `crs`, as anything `Proj.CRS` accepts,
e.g. `"EPSG:4326"` or a WKT string. With `crs=nothing` no `proj:` attributes are written.

## Reading a store

`Pyramid(path)` opens a multiscale store. The base and the order of the levels are taken from the
`derived_from` references and the scales of the `multiscales` layout, so stores written by other tools
can be read as well. Levels with a higher resolution than the base are skipped.
Stores written by earlier versions of PyramidScheme, with the base at the root of the store,
can still be read.

Reading from remote stores needs the corresponding storage package, e.g. `using Zarr` loads all of them.

## API

```@docs
buildpyramids
Base.write(::AbstractString, ::PyramidScheme.Pyramid)
PyramidScheme.MultiscalesStyle
ConventionStyle
CellCenterStyle
PyramidScheme.pyramid_factors
```
