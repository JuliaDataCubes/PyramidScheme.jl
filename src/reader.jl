# Reading pyramids from Zarr stores.
# Stores with a `multiscales` layout (https://github.com/zarr-conventions/multiscales) are read from
# the layout, older PyramidScheme stores keep the base at the root and the levels in groups "1", "2", ...

function _pyramid_zarr(path)
    g = zopen(path)
    multiscales = get(g.attrs, "multiscales", nothing)
    if multiscales isa AbstractDict && haskey(multiscales, "layout")
        read_multiscales(path, g, multiscales["layout"])
    else
        read_legacy_pyramid(path, g)
    end
end

function read_multiscales(path, g, layout)
    scales = layout_scales(layout)
    bases = [i for (i, l) in enumerate(layout) if !haskey(l, "derived_from")]
    isempty(bases) && throw(ArgumentError("The multiscales layout of $path has no level without derived_from"))
    finer = [l["asset"] for (l, s) in zip(layout, scales) if s < 1]
    isempty(finer) ||
        @warn "Skipping the levels $finer of $path, they have a higher resolution than the base"
    coarser = sort([i for i in eachindex(layout) if scales[i] > 1]; by=i -> scales[i])
    base = read_level(g, layout[first(bases)]["asset"])
    levels = [read_level(g, layout[i]["asset"]) for i in coarser]
    Pyramid(base, levels, Dict{String,Any}(g.attrs))
end

"""
    layout_scales(layout)

Return the factor by which each level of a multiscales `layout` is downsampled compared to the level
it is ultimately derived from, following the `derived_from` references and multiplying the spatial
`scale` of their `transform`s.
"""
function layout_scales(layout)
    byasset = Dict(l["asset"] => l for l in layout)
    function scale(l, depth)
        haskey(l, "derived_from") || return 1.0
        depth > length(layout) && throw(ArgumentError("The derived_from references of the multiscales layout form a cycle"))
        source = get(byasset, l["derived_from"], nothing)
        source === nothing &&
            throw(ArgumentError("Level $(l["asset"]) is derived from $(l["derived_from"]), which is not in the layout"))
        scale(source, depth + 1) * step_scale(l)
    end
    [scale(l, 0) for l in layout]
end

# The scale of the spatial dimensions, the other dimensions have a scale of one
function step_scale(l)
    s = get(get(l, "transform", Dict()), "scale", [1.0])
    i = findfirst(!=(1), s)
    i === nothing ? 1.0 : Float64(s[i])
end

function read_level(g, asset)
    node = zarr_node(g, asset)
    node isa ZarrCore.ZGroup && return Cube(open_dataset(node))
    parts = split(asset, '/')
    parentgroup = length(parts) == 1 ? g : zarr_node(g, join(parts[1:end-1], '/'))
    open_dataset(parentgroup).cubes[Symbol(last(parts))]
end

function zarr_node(g, asset)
    node = g
    for part in split(asset, '/')
        node = if node isa ZarrCore.ZGroup && haskey(node.groups, part)
            node.groups[part]
        elseif node isa ZarrCore.ZGroup && haskey(node.arrays, part)
            node.arrays[part]
        else
            throw(ArgumentError("The multiscales asset $asset is not in the store"))
        end
    end
    node
end

function read_legacy_pyramid(path, g)
    levelnames = filter(k -> occursin(r"^\d+$", k), collect(keys(g.groups)))
    base = Cube(path)
    clevels = [Cube(open_dataset(g[l])) for l in sort(levelnames; by=l -> parse(Int, l))]
    Pyramid(base, clevels, DD.metadata(base))
end

"""
    level_factors(pyramid; spatial_dims=SpatialDim)

Return the downsampling factor of every level of `pyramid`, starting with 1 for the base.
"""
function level_factors(pyramid; spatial_dims=SpatialDim)
    xbase = first(DD.dims(parent(pyramid), spatial_dims))
    [1; [level_factor(xbase, DD.dims(l, DD.basetypeof(xbase))) for l in pyramid.levels]]
end
