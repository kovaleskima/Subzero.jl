export BondSettings

"""
    BondSettings{FT}

Settings needed for bonds between floes.

## _Fields_ / _Keyword Arguments_
- `bonds_on::Bool`: if true, bonded floes are held together by springs (Default = false)
- `coastline_bonds::Bool`: if true, floes touching the coast or topography at the start get bonds to it (Default = false)
- `stiffness_ratio::FT`: bond stiffness is the contact stiffness divided by this ratio (Default = 25)
- `damping_ratio::FT`: fraction of critical damping applied to each bond (Default = 0.1)

Construct with `BondSettings([FT = Float64]; kwargs...)`.
"""
struct BondSettings{FT<:AbstractFloat}
    bonds_on::Bool
    coastline_bonds::Bool
    stiffness_ratio::FT
    damping_ratio::FT

    function BondSettings{FT}(
        bonds_on, coastline_bonds, stiffness_ratio, damping_ratio,
    ) where {FT<:AbstractFloat}
        if bonds_on
            if stiffness_ratio <= 0
                @warn "Bond stiffness ratio must be positive. Turning bonds off."
                bonds_on = false
            elseif damping_ratio < 0
                @warn "Bond damping ratio can't be negative. Turning bonds off."
                bonds_on = false
            end
        end
        if !bonds_on && coastline_bonds
            @warn "Coastline bonds can't be on without bonds. Turning coastline bonds off."
            coastline_bonds = false
        end
        new{FT}(bonds_on, coastline_bonds, stiffness_ratio, damping_ratio)
    end
end

function BondSettings(::Type{FT} = Float64;
    bonds_on = false, coastline_bonds = false,
    stiffness_ratio = 25, damping_ratio = 0.1,
) where {FT<:AbstractFloat}
    return BondSettings{FT}(bonds_on, coastline_bonds, FT(stiffness_ratio), FT(damping_ratio))
end

function Base.show(io::IO, bs::BondSettings{FT}) where {FT}
    print(io, "BondSettings{$FT}\n")
    for name in fieldnames(BondSettings)
        print(io, "  ⊢ ", name, " = ", getfield(bs, name), "\n")
    end
end