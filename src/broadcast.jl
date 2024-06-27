import Base.Broadcast: Broadcasted, BroadcastStyle

struct PyramidStyle <:Base.BroadcastStyle end
Base.BroadcastStyle(::Type{<:Pyramid}) = PyramidStyle()
Base.BroadcastStyle(::PyramidStyle, ::PyramidStyle) = PyramidStyle()

function Base.copy(bc::Broadcasted{PyramidStyle})
    bcf = Base.Broadcast.flatten(bc)
    inputs = bcf.args
    func = bcf.f
    numlevels = checklevelcompat(inputs)
    newlevels = map(0:numlevels) do l
        argslevel = levels.(inputs, (l,))
        argdata = getproperty.(argslevel, :data)
        newdata = func.(argdata...)
        newdimarr = DD.rebuild(first(argslevel), data=newdata, dims=DD.dims(first(argslevel)))
        newdimarr
    end
    #@show typeof(newlevels)
    base = newlevels[1]
    layers = newlevels[2:end]
    #@show typeof(base)
    #@show eltype(layers)
    Pyramid(base, layers, broadcast_metadata(inputs))
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