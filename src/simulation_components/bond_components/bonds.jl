# add BondSet struct, constructors, and alive flag
# the bond as seen from one of its two floes
struct FloeBond{FT<:AbstractFloat}
    partner_id::Int                  # stable `id` of the bonded floe
    r_self::SA.SVector{2,FT}         # bond center in THIS floe's body frame
    r_partner::SA.SVector{2,FT}      # same bond center in the partner's body frame
    L::FT                            # shared face length
end

# bond to the coast: fixed anchor, no partner floe
struct CoastBond{FT<:AbstractFloat}
    r_self::SA.SVector{2,FT}
    anchor::SA.SVector{2,FT}         # fixed world position
    L::FT
end

Base.@kwdef struct FloeBonds{FT<:AbstractFloat}
    floe::Vector{FloeBond{FT}} = FloeBond{FT}[]
    coast::Vector{CoastBond{FT}} = CoastBond{FT}[]
end

Base.isempty(b::FloeBonds) = isempty(b.floe) && isempty(b.coast)
nbonds(b::FloeBonds) = length(b.floe) + length(b.coast)