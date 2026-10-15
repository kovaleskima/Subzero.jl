export BondSettings

"""
    BondSettings{FT}

Settings needed for bonds between floes.

## _Fields_ / _Keyword Arguments_
- `bonds_on::Bool`: if true, bonded floes are held together by springs (Default = false)
- `coastline_bonds::Bool`: if true, bonds between floes and the coast/topography (added with
    `add_coast_bond!`) apply forces; if false, they are ignored and only floe-floe bonds act.
    Requires `bonds_on`. (Default = false)
- `stiffness_ratio::FT`: bond stiffness is the contact stiffness divided by this ratio (Default = 25)
- `damping_ratio::FT`: fraction of critical damping applied to each bond (Default = 0.1)
- `stress_nhistory::Int`: number of most recent timesteps of bond stress that are accumulated
    when checking bond fracture (Default = 10)
- `stress_decay::FT`: decay rate of the exponential weighting of accumulated bond stress. The
    stress from `k` timesteps ago has weight `exp(-stress_decay * k)`, so larger values put more
    weight on recent timesteps and `0` weights all `stress_nhistory` timesteps equally. Weights
    are normalized over the full `stress_nhistory` timesteps, so a newly created bond needs
    `stress_nhistory` timesteps of sustained stress to reach its full accumulated stress.
    (Default = 0.2, so the oldest timestep has ~1/6 the weight of the newest)

Construct with `BondSettings([FT = Float64]; kwargs...)`.
"""
struct BondSettings{FT<:AbstractFloat}
    bonds_on::Bool
    coastline_bonds::Bool
    stiffness_ratio::FT
    damping_ratio::FT
    stress_nhistory::Int
    stress_decay::FT

    function BondSettings{FT}(
        bonds_on, coastline_bonds, stiffness_ratio, damping_ratio,
        stress_nhistory = 10, stress_decay = 0.2,
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
        if stress_nhistory < 1
            @warn "Bond stress must be accumulated over at least one timestep. Resetting \
                stress_nhistory to default value of 10."
            stress_nhistory = 10
        end
        if stress_decay < 0
            @warn "Bond stress decay can't be negative. Resetting stress_decay to default \
                value of 0.2."
            stress_decay = 0.2
        end
        if !bonds_on && coastline_bonds
            @warn "Coastline bonds can't be on without bonds. Turning coastline bonds off."
            coastline_bonds = false
        end
        new{FT}(bonds_on, coastline_bonds, stiffness_ratio, damping_ratio,
            stress_nhistory, stress_decay)
    end
end

function BondSettings(::Type{FT} = Float64;
    bonds_on = false, coastline_bonds = false,
    stiffness_ratio = 25, damping_ratio = 0.1, stress_nhistory = 10, stress_decay = 0.2,
) where {FT<:AbstractFloat}
    return BondSettings{FT}(bonds_on, coastline_bonds, FT(stiffness_ratio), FT(damping_ratio),
        stress_nhistory, FT(stress_decay))
end

function Base.show(io::IO, bs::BondSettings{FT}) where {FT}
    print(io, "BondSettings{$FT}\n")
    for name in fieldnames(BondSettings)
        print(io, "  ⊢ ", name, " = ", getfield(bs, name), "\n")
    end
end