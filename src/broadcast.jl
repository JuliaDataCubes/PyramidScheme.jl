import Base.Broadcast: Broadcasted, BroadcastStyle, DefaultArrayStyle, AbstractArrayStyle, Style

struct PyramidStyle{S <: BroadcastStyle} <: AbstractArrayStyle{Any} end
PyramidStyle(::S) where {S} = PyramidStyle{S}()
PyramidStyle(::S, ::Val{N}) where {S,N} = PyramidStyle(S(Val(N)))
PyramidStyle(::Val{N}) where N = PyramidStyle{DefaultArrayStyle{N}}()
function BroadcastStyle(::Type{<:Pyramid{T,N,D,A, B, L, Me}}) where {T,N,D,A, B, L, Me}
    inner_style = typeof(BroadcastStyle(B))
    return PyramidStyle{inner_style}()
end

function PyramidStyle(a::BroadcastStyle, b::BroadcastStyle)
    inner_style = BroadcastStyle(a, b)
    # if the inner style is Unknown then so is the outer style
    if inner_style isa Unknown
        return Unknown()
    else
        return PyramidStyle(inner_style)
    end
end


BroadcastStyle(::PyramidStyle, ::Base.Broadcast.Unknown) = Unknown()
BroadcastStyle(::Base.Broadcast.Unknown, ::PyramidStyle) = Unknown()
BroadcastStyle(::PyramidStyle{A}, ::PyramidStyle{B}) where {A, B} = PyramidStyle(A(), B())
BroadcastStyle(::PyramidStyle{A}, b::Style) where {A} = PyramidStyle(A(), b)
BroadcastStyle(a::Style, ::PyramidStyle{B}) where {B} = PyramidStyle(a, B())
BroadcastStyle(::PyramidStyle{A}, b::Style{Tuple}) where {A} = PyramidStyle(A(), b)
BroadcastStyle(a::Style{Tuple}, ::PyramidStyle{B}) where {B} = PyramidStyle(a, B())
BroadcastStyle(a::PyramidStyle, ::DD.DimensionalStyle) = a
BroadcastStyle(a::PyramidStyle, ::DiskArrays.ChunkStyle) = a

function Base.copy(bc::Broadcasted{<:PyramidStyle})
    bcf = Base.Broadcast.flatten(bc)
    args = bcf.args
    func = bcf.f
    arrs = Flatten.flatten(args, AbstractArray, Pyramid)
    isempty(arrs) || 
        throw(ArgumentError("Cannot broadcast a Pyramid with a regular array. Convert your input array to a Pyramid using `Pyramid(array)`"))
    pyrs = Flatten.flatten(args, Pyramid, Array)
    numlevels = checklevelcompat(pyrs)
    newlevels = map(0:numlevels) do l
        pyrslevel = levels.(pyrs, (l,))
        argslevel = Flatten.reconstruct(args, pyrslevel, Pyramid, Array) 
        func.(argslevel...)
    end
    #@show typeof(newlevels)
    base = newlevels[1]
    layers = newlevels[2:end]
    Pyramid(base, layers, broadcast_metadata(args))
end

#function Base.copyto!(pyr::Pyramid, bc::Broadcasted{PyramidStyle}) = pyr

# TODO Implement proper metadata handling for broadcast
broadcast_metadata(inputs) = Dict()

function checklevelcompat(inputs)
    levfirstpyr = nlevels(first(inputs))
    for inp in inputs
        if nlevels(inp) != levfirstpyr
            error("Pyramids should have the same number of levels got $levfirstpyr and $(nlevels(inp))")
        end
    end
    return levfirstpyr
end