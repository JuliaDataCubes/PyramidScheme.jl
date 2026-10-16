# Grids of the pyramid levels.
# A level downsampled by `factor` keeps the origin (top-left corner) of the base and has cells
# `factor` times larger, see `spatial_transform` in conventions.jl.

"""
    pyramid_factors(; nlevels=nothing, factors=nothing, default_nlevels=0)

Return the cumulative downsampling factors of all levels, starting with `1` for the base.
Either give the number of aggregated levels `nlevels`, which are downsampled by powers of two,
or the `factors` directly, e.g. `[1, 2, 8]`. Without either, `default_nlevels` levels are used.
"""
function pyramid_factors(; nlevels=nothing, factors=nothing, default_nlevels=0)
    nlevels === nothing || factors === nothing ||
        throw(ArgumentError("Give either nlevels or factors, not both"))
    factors === nothing || return check_factors(collect(Int, factors))
    n = nlevels === nothing ? default_nlevels : nlevels
    n >= 0 || throw(ArgumentError("nlevels must not be negative, got $n"))
    [2^i for i in 0:n]
end

"""
    level_length(n, factor, edges)

Length of a dimension of length `n` after downsampling by `factor`.
With `edges=:keep` a partial window at the end becomes its own cell,
with `edges=:drop` it is left out.
"""
function level_length(n, factor, edges)
    len = if edges === :keep
        cld(n, factor)
    elseif edges === :drop
        n ÷ factor
    else
        throw(ArgumentError("edges must be :keep or :drop, got $(repr(edges))"))
    end
    len >= 1 || throw(ArgumentError("A dimension of length $n can not be downsampled by $factor with edges=:drop"))
    len
end

"""
    coarsen_dim(d, factor; edges=:keep)

Return the dimension `d` downsampled by `factor`.
The new cells are `factor` times larger and start at the same corner as the cells of `d`,
the lookup keeps its type, order and locus.
"""
function coarsen_dim(d, factor; edges=:keep)
    n = level_length(length(d), factor, edges)
    factor == 1 && return d
    res = cellsize(d)
    newres = res * factor
    corner = first_cell_center(d, res) - res / 2
    centers = range(corner + newres / 2; step=newres, length=n)
    l = DD.lookup(d)
    l isa DD.AbstractSampled || return DD.rebuild(d, centers)
    values = if DD.sampling(l) isa DD.Intervals
        locus = DD.locus(l)
        locus isa DD.Start ? centers .- newres / 2 : locus isa DD.End ? centers .+ newres / 2 : centers
    else
        centers
    end
    DD.rebuild(d, DD.rebuild(l; data=values, span=DD.Regular(newres)))
end

"""
    level_dims(dims, factor; edges=:keep, spatial_dims=SpatialDim)

Return `dims` with the spatial dimensions downsampled by `factor`, see [`coarsen_dim`](@ref).
Other dimensions are returned unchanged.
"""
level_dims(dims, factor; edges=:keep, spatial_dims=SpatialDim) =
    map(d -> isspatial(dims, d, spatial_dims) ? coarsen_dim(d, factor; edges) : d, dims)

isspatial(dims, d, spatial_dims) = DD.hasdim(DD.dims(dims, spatial_dims), d)

"""
    level_factor(basedim, leveldim)

Return the factor by which `leveldim` is downsampled compared to `basedim`.
"""
function level_factor(basedim, leveldim)
    nbase, nlevel = length(basedim), length(leveldim)
    nbase == nlevel && return 1
    nlevel == 1 && return nbase
    resbase, reslevel = regular_step(basedim), regular_step(leveldim)
    resbase === nothing || reslevel === nothing ? round(Int, nbase / nlevel) : round(Int, reslevel / resbase)
end

"""
    levelindex(factor, i, n)

Return the index into a level of length `n` that is downsampled by `factor`,
for the index `i` into the base.
"""
levelindex(factor, i::Integer, n) = min((i - 1) ÷ factor + 1, n)
levelindex(factor, i::AbstractUnitRange, n) = levelindex(factor, first(i), n):levelindex(factor, last(i), n)
levelindex(factor, ::Colon, n) = Colon()
