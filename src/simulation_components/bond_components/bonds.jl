# add BondSet struct, constructors, and alive flag

struct BondSet{FT<:AbstractFloat}
    i::Vector{Int}                    # stable floe IDs
    j::Vector{Int}
    ri::Vector{SA.SVector{2,FT}}         # bond center, floe i body frame
    rj::Vector{SA.SVector{2,FT}}
    L::Vector{FT}                     # shared face length
    alive::Vector{Bool}
end

Base.@kwdef struct BondParameters{FT<:AbstractFloat}
    stiffness_ratio::FT = 25     # k_bond = k_series / stiffness_ratio (MATLAB's "/25")
    damping_ratio::FT   = 0.1    # ζ, fraction of critical damping (arbitrary for now)
    coastline_bonds::Bool = false
end