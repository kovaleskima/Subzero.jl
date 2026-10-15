# add BondSet struct, constructors, and alive flag
# the bond as seen from one of its two floes
struct FloeBond{FT<:AbstractFloat}
    partner_id::Int
    attach_self::SA.SVector{2,FT}      # vector from floe centroid to bond anchor on the floe face 
    attach_partner::SA.SVector{2,FT}   # same thing but for partner floe
    normal_self::SA.SVector{2,FT}      # this floe's unit outward face normal, body frame
    L::FT                              # length of shared face
    stress_history::Vector{SA.SVector{2,FT}}  # (σn, τ) of recent timesteps, newest first
end

FloeBond{FT}(partner_id, attach_self, attach_partner, normal_self, L) where {FT} =
    FloeBond{FT}(partner_id, attach_self, attach_partner, normal_self, L, SA.SVector{2,FT}[])

# bond to the coast: fixed anchor, no partner floe
struct CoastBond{FT<:AbstractFloat}
    attach_self::SA.SVector{2,FT}      # same as above
    anchor::SA.SVector{2,FT}           # fixed world position
    normal_self::SA.SVector{2,FT}
    L::FT
    stress_history::Vector{SA.SVector{2,FT}}  # same as above
end

CoastBond{FT}(attach_self, anchor, normal_self, L) where {FT} =
    CoastBond{FT}(attach_self, anchor, normal_self, L, SA.SVector{2,FT}[])

Base.@kwdef mutable struct FloeBonds{FT<:AbstractFloat}
    floe::Vector{FloeBond{FT}} = FloeBond{FT}[]
    coast::Vector{CoastBond{FT}} = CoastBond{FT}[]
    fx::FT = zero(FT)
    fy::FT = zero(FT)
    trq::FT = zero(FT)
end

Base.isempty(b::FloeBonds) = isempty(b.floe) && isempty(b.coast)
nbonds(b::FloeBonds) = length(b.floe) + length(b.coast)