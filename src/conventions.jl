# Metadata for the zarr-conventions used by multiscale stores:
# https://github.com/zarr-conventions/multiscales
# https://github.com/zarr-conventions/proj
# https://github.com/zarr-conventions/spatial

"""
    Convention(repo, tag, uuid, name, description)

A Zarr convention pinned to one tag of its repository in the zarr-conventions organisation.
The schema and spec URLs are derived from `repo` and `tag`.
"""
struct Convention
    repo::String
    tag::String
    uuid::String
    name::String
    description::String
end

schema_url(c::Convention) =
    "https://raw.githubusercontent.com/zarr-conventions/$(c.repo)/refs/tags/$(c.tag)/schema.json"
spec_url(c::Convention) = "https://github.com/zarr-conventions/$(c.repo)/blob/$(c.tag)/README.md"

"""
    convention_attrs(c::Convention)

Return the entry for `c` in the `zarr_conventions` attribute.
"""
convention_attrs(c::Convention) = Dict{String,Any}(
    "schema_url" => schema_url(c),
    "spec_url" => spec_url(c),
    "uuid" => c.uuid,
    "name" => c.name,
    "description" => c.description,
)

const MULTISCALES_CONVENTION = Convention("multiscales", "v0.1",
    "d35379db-88df-4056-af3a-620245f8e347", "multiscales", "Multiscale layout of zarr datasets")
const PROJ_CONVENTION = Convention("proj", "v0.1",
    "f17cb550-5864-4468-aeb7-f3180cfb622f", "proj",
    "Coordinate reference system information for geospatial data")
const SPATIAL_CONVENTION = Convention("spatial", "v0.1",
    "689b58e2-cf7b-45e0-9fff-9cfc0883d6b4", "spatial", "Spatial coordinate information")

"""
    MultiscalesStyle

How the `multiscales` metadata is written where existing writers differ in their reading
of the convention:

- [`ConventionStyle`](@ref) follows the README and examples of the multiscales convention.
- [`CellCenterStyle`](@ref) measures translations between cell centers, as topozarr does,
  for clients built against that output.
"""
abstract type MultiscalesStyle end

"""
    ConventionStyle()

Write the `multiscales` metadata as in the README and examples of the multiscales convention:
the `translation` between levels is zero, level 0 has an identity `transform`,
and resampling methods use the names listed by the convention (`"average"` for `mean`).
"""
struct ConventionStyle <: MultiscalesStyle end

"""
    CellCenterStyle()

Write the `multiscales` metadata with translations measured between cell centers:
the `translation` between levels is `(step - 1) / 2` along the spatial dimensions,
level 0 has no `transform`, and `mean` is named `"mean"`.
This is the metadata written by topozarr (https://github.com/carbonplan/topozarr).
"""
struct CellCenterStyle <: MultiscalesStyle end

level_translation(::ConventionStyle, step) = 0.0
level_translation(::CellCenterStyle, step) = (step - 1) / 2

has_base_transform(::ConventionStyle) = true
has_base_transform(::CellCenterStyle) = false

const CELLCENTER_METHOD_NAMES = Dict("average" => "mean", "med" => "median")

resampling_method_name(::ConventionStyle, f) = resampling_method_name(f)
function resampling_method_name(::CellCenterStyle, f)
    name = resampling_method_name(f)
    get(CELLCENTER_METHOD_NAMES, name, name)
end

"""
    resampling_method_name([style,] f)

Return the name of the aggregation function `f` as used in the `resampling_method` field of the
multiscales convention, e.g. `"average"` for `mean`.
`mean ∘ skipmissing` is named after `mean`. Unknown functions fall back to `string(f)`.
With a [`MultiscalesStyle`](@ref) the name follows that style.
"""
resampling_method_name(f) = string(f)
resampling_method_name(f::ComposedFunction) =
    f.inner === skipmissing ? resampling_method_name(f.outer) : string(f)
resampling_method_name(::typeof(mean)) = "average"
resampling_method_name(::typeof(median)) = "med"
resampling_method_name(::typeof(sum)) = "sum"
resampling_method_name(::typeof(maximum)) = "max"
resampling_method_name(::typeof(minimum)) = "min"
resampling_method_name(::typeof(first)) = "nearest"
resampling_method_name(::ESALCMode) = "mode"

"""
    proj_attrs(crs)

Return the `proj:` attributes for `crs`, which can be anything `Proj.CRS` accepts,
e.g. `"EPSG:4326"`, a WKT string or a GeoFormatTypes CRS.
`proj:wkt2` is always written, on a single line. `proj:code` is only written when the CRS has an
authority code, or one is identified with full confidence.
"""
function proj_attrs(crs)
    pcrs = crs isa Proj.CRS ? crs : Proj.CRS(crs)
    attrs = Dict{String,Any}("proj:wkt2" => wkt2(pcrs))
    code = authority_code(pcrs)
    code === nothing || (attrs["proj:code"] = code)
    attrs
end

function wkt2(pcrs::Proj.CRS)
    option = "MULTILINE=NO"
    GC.@preserve option begin
        # NULL terminated list of PROJ options
        options = [Base.unsafe_convert(Cstring, option), reinterpret(Cstring, C_NULL)]
        GC.@preserve options Proj.proj_as_wkt(pcrs, Proj.PJ_WKT2_2019, pointer(options))
    end
end

function authority_code(pcrs::Proj.CRS)
    auth, code = Proj.proj_get_id_auth_name(pcrs), Proj.proj_get_id_code(pcrs)
    if auth === nothing || code === nothing
        matches = filter(m -> m.confidence == 100, Proj.identify(pcrs))
        length(matches) == 1 || return nothing
        auth, code = Proj.proj_get_id_auth_name(only(matches).crs), Proj.proj_get_id_code(only(matches).crs)
        (auth === nothing || code === nothing) && return nothing
    end
    "$auth:$code"
end

"""
    cellsize(d; fallback=nothing)

Return the signed distance between neighbouring cells of the dimension `d`.
Throws an `ArgumentError` when `d` is not regularly spaced.
A dimension of length one has no spacing, so `fallback` is returned instead.
"""
function cellsize(d; fallback=nothing)
    if length(d) < 2
        fallback === nothing &&
            throw(ArgumentError("Can not infer the cell size of dimension $(DD.name(d)) from a single value"))
        return Float64(fallback)
    end
    step = regular_step(d)
    step === nothing && throw(ArgumentError("Dimension $(DD.name(d)) is not regularly spaced"))
    step
end

"""
    regular_step(d)

Return the distance between neighbouring values of the dimension `d`,
or `nothing` when `d` is shorter than two or not regularly spaced.
"""
function regular_step(d)
    length(d) < 2 && return nothing
    steps = diff(Float64.(collect(DD.lookup(d))))
    all(s -> isapprox(s, first(steps); rtol=1e-5), steps) ? first(steps) : nothing
end

"""
    first_cell_center(d, res)

Return the coordinate of the center of the first cell of `d`, with cell size `res`.
Point lookups are taken as cell centers, interval lookups are shifted according to their locus.
"""
function first_cell_center(d, res)
    l = DD.lookup(d)
    v = Float64(first(l))
    l isa DD.AbstractSampled && DD.sampling(l) isa DD.Intervals || return v
    locus = DD.locus(l)
    locus isa DD.Start ? v + res / 2 : locus isa DD.End ? v - res / 2 : v
end

"""
    spatial_transform(xdim, ydim; fallback_res=nothing)

Return the affine transform `[a, b, c, d, e, f]` of the spatial convention for the grid given by
`xdim` and `ydim`. Index `(0, 0)` is the top-left corner of the first cell.
`fallback_res` is an `(xres, yres)` tuple used for dimensions of length one.
"""
function spatial_transform(xdim, ydim; fallback_res=nothing)
    xres = cellsize(xdim; fallback=fallback_res === nothing ? nothing : fallback_res[1])
    yres = cellsize(ydim; fallback=fallback_res === nothing ? nothing : fallback_res[2])
    c = first_cell_center(xdim, xres) - xres / 2
    f = first_cell_center(ydim, yres) - yres / 2
    [xres, 0.0, c, 0.0, yres, f]
end

"""
    spatial_bbox(transform, shape)

Return the bounding box `[xmin, ymin, xmax, ymax]` of a grid with affine `transform`
and spatial `shape` `[height, width]`.
"""
function spatial_bbox(transform, shape)
    xres, _, c, _, yres, f = transform
    ny, nx = shape
    xs, ys = (c, c + xres * nx), (f, f + yres * ny)
    [minimum(xs), minimum(ys), maximum(xs), maximum(ys)]
end

"""
    spatial_attrs(xdim, ydim)

Return the `spatial:` attributes describing the grid given by `xdim` and `ydim`.
"""
function spatial_attrs(xdim, ydim)
    transform = spatial_transform(xdim, ydim)
    shape = [length(ydim), length(xdim)]
    Dict{String,Any}(
        "spatial:dimensions" => [string(DD.name(ydim)), string(DD.name(xdim))],
        "spatial:registration" => "pixel",
        "spatial:transform" => transform,
        "spatial:bbox" => spatial_bbox(transform, shape),
        "spatial:shape" => shape,
    )
end

"""
    check_factors(factors)

Check that the cumulative downsampling `factors` start at one, increase strictly
and that every factor is a multiple of the previous one.
"""
function check_factors(factors)
    isempty(factors) && throw(ArgumentError("factors must not be empty"))
    first(factors) == 1 || throw(ArgumentError("factors must start at 1, got $factors"))
    for (prev, cur) in zip(factors[1:end-1], factors[2:end])
        cur > prev || throw(ArgumentError("factors must increase strictly, got $factors"))
        cur % prev == 0 ||
            throw(ArgumentError("each factor must be a multiple of the previous one, got $factors"))
    end
    factors
end

"""
    multiscales_attrs(leveldims, factors; resampling_method, recursive=true, spatial_dims=SpatialDim, style=ConventionStyle())

Return the value of the `multiscales` attribute for a pyramid whose level `i` (starting at 0)
has the dimensions `leveldims[i+1]` and is downsampled by the cumulative factor `factors[i+1]`
along the spatial dimensions.

Every level is stored in the child group named after its index. With `recursive=true` each level
is derived from the previous one, otherwise all levels are derived from level 0.
`scale` and `translation` list every dimension in stored (row-major) order,
the reverse of the Julia order.
`style` selects how translations, the `transform` of level 0 and the resampling method name are
written, see [`MultiscalesStyle`](@ref).
"""
function multiscales_attrs(leveldims, factors; resampling_method, recursive=true, spatial_dims=SpatialDim,
    style::MultiscalesStyle=ConventionStyle())
    length(leveldims) == length(factors) ||
        throw(ArgumentError("Got $(length(leveldims)) levels but $(length(factors)) factors"))
    check_factors(factors)
    method = resampling_method_name(style, resampling_method)
    xdim0, ydim0 = DD.dims(first(leveldims), spatial_dims)
    res0 = cellsize(xdim0), cellsize(ydim0)
    layout = map(enumerate(zip(leveldims, factors))) do (i, (dims, factor))
        source = i == 1 ? nothing : recursive ? i - 1 : 1
        step = source === nothing ? 1 : factor ÷ factors[source]
        isspatial = map(d -> DD.hasdim(DD.dims(dims, spatial_dims), d), reverse(dims))
        xdim, ydim = DD.dims(dims, spatial_dims)
        level = Dict{String,Any}(
            "asset" => string(i - 1),
            "spatial:transform" => spatial_transform(xdim, ydim; fallback_res=res0 .* factor),
            "spatial:shape" => [length(ydim), length(xdim)],
        )
        if source !== nothing || has_base_transform(style)
            level["transform"] = Dict{String,Any}(
                "scale" => [s ? Float64(step) : 1.0 for s in isspatial],
                "translation" => [s ? level_translation(style, step) : 0.0 for s in isspatial],
            )
        end
        if source !== nothing
            level["derived_from"] = string(source - 1)
            level["resampling_method"] = method
        end
        level
    end
    Dict{String,Any}("layout" => collect(layout), "resampling_method" => method)
end

"""
    multiscales_group_attrs(leveldims, factors; crs=nothing, resampling_method, recursive=true, spatial_dims=SpatialDim, style=ConventionStyle())

Return the attributes of the root group of a multiscale store: the `zarr_conventions` registry,
the `multiscales` layout, the `spatial:` attributes of level 0 and,
when `crs` is given, the `proj:` attributes.
See [`multiscales_attrs`](@ref) for the meaning of the other arguments.
"""
function multiscales_group_attrs(leveldims, factors; crs=nothing, resampling_method, recursive=true,
    spatial_dims=SpatialDim, style::MultiscalesStyle=ConventionStyle())
    conventions = crs === nothing ? [MULTISCALES_CONVENTION, SPATIAL_CONVENTION] :
                  [MULTISCALES_CONVENTION, PROJ_CONVENTION, SPATIAL_CONVENTION]
    attrs = Dict{String,Any}("zarr_conventions" => convention_attrs.(conventions))
    crs === nothing || merge!(attrs, proj_attrs(crs))
    merge!(attrs, spatial_attrs(DD.dims(first(leveldims), spatial_dims)...))
    attrs["multiscales"] = multiscales_attrs(leveldims, factors; resampling_method, recursive, spatial_dims, style)
    attrs
end
