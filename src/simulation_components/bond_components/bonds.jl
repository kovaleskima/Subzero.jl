# add BondSet struct, constructors, and alive flag

struct BondSet{FT<:AbstractFloat}
    i::Vector{Int}                    # stable floe IDs
    j::Vector{Int}
    ri::Vector{SVector{2,FT}}         # bond center, floe i body frame
    rj::Vector{SVector{2,FT}}
    L::Vector{FT}                     # shared face length
    alive::Vector{Bool}
end