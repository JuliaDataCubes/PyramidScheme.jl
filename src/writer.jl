# Writing pyramids as Zarr v3 stores following the multiscales convention:
# the root group only holds the metadata, every level including the base is a child group
# named after its index.

default_compressor() = ZarrZstd.ZstdCompressor()

"""
    buildpyramids(src, dest; nlevels=nothing, factors=nothing, resampling_method=mean ∘ skipmissing,
                  recursive=true, edges=:keep, crs=:auto, spatial_dims=SpatialDim, tilesize=256, chunksize=1024,
                  compressor=ZstdCompressor(), consolidated=true, overwrite=false, runner=LocalRunner)

Build a multiscale Zarr v3 store at `dest` from `src` and return it as a [`Pyramid`](@ref).
`src` is the path of a dataset that YAXArrays can open, or an `AbstractDimArray`. It is not modified.

The store follows the multiscales, spatial and proj conventions of
https://github.com/zarr-conventions: the root group holds the metadata, `src` is copied to the child
group `0` and the aggregated levels are written to the groups `1`, `2`, ...

# Keywords
- `nlevels`, `factors`: the levels to build, see [`pyramid_factors`](@ref). Without either, enough levels
  are built for the coarsest one to fit into `tilesize`.
- `resampling_method`: function aggregating a window of cells, e.g. `mean`, `sum` or `mean ∘ skipmissing`.
  The default `mean ∘ skipmissing` leaves missing values out, a window of only missing values gives `NaN`.
- `recursive`: compute every level from the previous one instead of from `src`.
- `edges`: `:keep` turns partial windows at the end of a dimension into a cell, `:drop` leaves them out.
- `crs`: the coordinate reference system written to the `proj:` attributes, anything `Proj.CRS` accepts.
  `:auto` takes it from the lookups of the spatial dimensions if they have one,
  `nothing` writes no `proj:` attributes.
- `spatial_dims`: the X and Y dimensions, in this order.
- `chunksize`: chunk length along the spatial dimensions. Other dimensions use chunks of length one.
- `compressor`: compressor of the data arrays.
- `consolidated`: write consolidated metadata into the root group.
- `overwrite`: replace `dest` if it exists.
"""
buildpyramids(src::AbstractString, dest::AbstractString; kwargs...) = buildpyramids(Cube(src), dest; kwargs...)

function buildpyramids(src::DD.AbstractDimArray, dest::AbstractString; nlevels=nothing, factors=nothing,
    resampling_method=mean ∘ skipmissing, recursive=true, edges=:keep, crs=:auto, spatial_dims=SpatialDim,
    tilesize=256, chunksize=1024, compressor=default_compressor(), consolidated=true, overwrite=false, runner=LocalRunner)
    xy = DD.dims(src, spatial_dims)
    length(xy) == 2 || throw(ArgumentError("Expected two spatial dimensions got $xy"))
    factors = pyramid_factors(; nlevels, factors, default_nlevels=compute_nlevels(src, tilesize; spatial_dims))
    leveldims = [level_dims(DD.dims(src), f; edges, spatial_dims) for f in factors]
    crs = crs === :auto ? default_crs(src, spatial_dims) : crs
    attrs = multiscales_group_attrs(leveldims, factors; crs, resampling_method, recursive, spatial_dims)
    layername = layer_name(src)
    levelattrs = level_group_attrs(attrs)

    create_multiscales_group(dest, attrs; overwrite)
    write_level(src, levelpath(dest, 0), spatial_dims; layername, chunksize, compressor, attrs=levelattrs)
    if length(factors) > 1
        # The levels are computed with partial windows and trimmed for edges=:drop,
        # the trimmed cells do not contribute to coarser levels.
        compute_levels(src, factors, resampling_method; recursive, spatial_dims, runner) do i, d, t
            output_zarr(d, t, levelpath(dest, i), spatial_dims; layername, chunksize, compressor, attrs=levelattrs)
        end
        if edges === :drop
            for i in 1:length(factors)-1
                trim_level!(levelpath(dest, i), leveldims[i+1])
            end
        end
    end
    consolidated && ZarrCore.consolidate_metadata(zopen(dest, "w"))
    open_written_pyramid(dest, length(factors), attrs)
end

buildpyramids(path::AbstractString; kwargs...) = throw(ArgumentError("""
    Building a pyramid inside the source dataset is no longer supported, because the multiscales
    convention stores every level, including the base, in its own child group.
    Use `buildpyramids(src, dest)` to write the pyramid to a new store at `dest`."""))

"""
    write(dest, pyramid::Pyramid; resampling_method, recursive=true, crs=:auto,
          spatial_dims=SpatialDim, chunksize=1024, compressor=ZstdCompressor(), consolidated=true, overwrite=false)

Write `pyramid` to a multiscale Zarr v3 store at `dest`, see [`buildpyramids`](@ref) for the layout
and the keywords. `resampling_method` defaults to the method recorded in the metadata of `pyramid`.
"""
function Base.write(dest::AbstractString, pyramid::Pyramid; resampling_method=recorded_method(pyramid),
    recursive=true, crs=:auto, spatial_dims=SpatialDim, chunksize=1024, compressor=default_compressor(),
    consolidated=true, overwrite=false)
    base = parent(pyramid)
    alllevels = [base, pyramid.levels...]
    factors = level_factors(pyramid; spatial_dims)
    crs = crs === :auto ? default_crs(base, spatial_dims) : crs
    attrs = multiscales_group_attrs(DD.dims.(alllevels), factors; crs, resampling_method, recursive, spatial_dims)
    layername = layer_name(base)
    levelattrs = level_group_attrs(attrs)

    create_multiscales_group(dest, attrs; overwrite)
    for (i, level) in enumerate(alllevels)
        write_level(level, levelpath(dest, i - 1), spatial_dims; layername, chunksize, compressor, attrs=levelattrs)
    end
    consolidated && ZarrCore.consolidate_metadata(zopen(dest, "w"))
    open_written_pyramid(dest, length(alllevels), attrs)
end

function recorded_method(pyramid)
    md = DD.metadata(pyramid)
    md isa AbstractDict && haskey(md, "resampling_method") && return md["resampling_method"]
    throw(ArgumentError("The pyramid does not record its resampling method, pass it as `resampling_method`"))
end

levelpath(dest, i) = joinpath(dest, string(i))

"""
    default_crs(a, spatial_dims)

Return the CRS stored in the lookups of the spatial dimensions of `a`, e.g. by Rasters.jl,
or `nothing`.
"""
function default_crs(a, spatial_dims)
    for d in DD.dims(a, spatial_dims)
        l = DD.lookup(d)
        crs = hasproperty(l, :crs) ? getproperty(l, :crs) : nothing
        crs === nothing || return crs
    end
    nothing
end

function layer_name(a)
    n = DD.name(a)
    n isa Symbol && !isempty(string(n)) ? string(n) : "layer"
end

function create_multiscales_group(dest, attrs; overwrite=false)
    if ispath(dest)
        overwrite || throw(ArgumentError("$dest already exists, set `overwrite=true` to replace it"))
        rm(dest; recursive=true)
    end
    zgroup(ZarrCore.storefromstring(dest, true)..., 3; attrs)
end

level_chunks(dims, spatial_dims, chunksize) =
    map(d -> isspatial(dims, d, spatial_dims) ? min(chunksize, length(d)) : 1, dims)

"""
    write_level(a, path, spatial_dims; layername, chunksize, compressor, attrs)

Write the array `a` as a Zarr v3 dataset to `path`, with the group attributes `attrs`.
"""
function write_level(a, path, spatial_dims; layername, chunksize, compressor, attrs)
    yax = a isa YAXArray ? a : YAXArray(DD.dims(a), parent(a))
    chunked = setchunks(yax, level_chunks(DD.dims(a), spatial_dims, chunksize))
    ds = to_dataset(chunked; layername)
    merge!(ds.properties, attrs)
    savedataset(ds; path, driver=:zarr, zarr_format=3, compressor)
end

"""
    output_zarr(leveldims, t, path, spatial_dims; layername, chunksize, compressor, attrs)

Construct an empty Zarr v3 dataset at `path` for a pyramid level with the dimensions `leveldims`,
element type `t` and group attributes `attrs`, and return the array of the level.
"""
function output_zarr(leveldims, t, path, spatial_dims; layername, chunksize, compressor, attrs)
    yax = YAXArray(leveldims, Zeros(t, length.(leveldims)...))
    chunked = setchunks(yax, level_chunks(leveldims, spatial_dims, chunksize))
    ds = to_dataset(chunked; layername)
    merge!(ds.properties, attrs)
    dssaved = savedataset(ds; path, skeleton=true, driver=:zarr, zarr_format=3, compressor)
    dssaved.cubes[Symbol(layername)].data
end

"""
    trim_level!(path, leveldims)

Shrink the arrays of the level at `path` to the lengths of `leveldims`.
"""
function trim_level!(path, leveldims)
    lengths = Dict(string(DD.name(d)) => length(d) for d in leveldims)
    g = zopen(path, "w")
    for a in values(g.arrays)
        newsize = map((name, n) -> get(lengths, name, n), ZarrCore.dimension_names(a), size(a))
        newsize == size(a) || resize!(a, newsize)
    end
end

function open_written_pyramid(dest, nlevels, attrs)
    base = Cube(levelpath(dest, 0))
    levels = [Cube(levelpath(dest, i)) for i in 1:nlevels-1]
    Pyramid(base, levels, attrs)
end
