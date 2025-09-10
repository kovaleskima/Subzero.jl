export AbstractBond, NoBond, Bond

abstract type AbstractBond end

struct NoBond <: AbstractBond end

mutable struct Bond{FT<:AbstractFloat}<:AbstractBond
    floe_id::Int
    neighbor_id::Int
end

"""
    Bond(::Type{FT}, args...)

A float type FT can be provided as the first argument of any Bond
constructor. A Bond of type FT will be created by passing all
other arguments to the correct constructor. 
"""
Bond(::Type{FT}, args...) where {FT <: AbstractFloat}=
    Bond{FT}(args...)

"""
    Bond(args...)

If a type isn't specified, Bond will be of type Float64 and the
correct constructor will be called with all other arguments.
"""
Bond(args...) = Bond{Float64}(args...)