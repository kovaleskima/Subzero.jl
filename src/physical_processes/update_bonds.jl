function timestep_bonds!(floes::StructArray{<:Floe{FT}}, domain, consts, bs, Δt) where {FT}
    idx_of = Dict{Int,Int}(floes.id[i] => i for i in eachindex(floes))
    Lx, Ly = _periodic_lengths(domain, FT)       # Inf on non-periodic axes
    Threads.@threads for i in eachindex(floes)
        bonds = floes.bonds[i]
        bonds.fx = bonds.fy = bonds.trq = zero(FT)     # always reset
        isempty(bonds) && continue
        xi = SA.SVector{2,FT}(floes.centroid[i][1], floes.centroid[i][2])
        Ui = SA.SVector{2,FT}(floes.u[i], floes.v[i])
        for b in bonds.floe
            j = get(idx_of, b.partner_id, 0)
            j == 0 && continue                          # dangling; pruned elsewhere
            xj0 = SA.SVector{2,FT}(floes.centroid[j][1], floes.centroid[j][2])
            xj = xi + _wrap(xj0 - xi, Lx, Ly)           # minimum image
            Uj = SA.SVector{2,FT}(floes.u[j], floes.v[j])
            ih, jh = floes.height[i], floes.height[j]
            ir, jr = sqrt(floes.area[i]), sqrt(floes.area[j])
            k = consts.E * (ih*jh) / (ih*jr + jh*ir) / bs.stiffness_ratio
            F, τ, _ = bond_force(xi, floes.α[i], b.attach_self, xj, floes.α[j], b.attach_partner, k, b.L)
            mi, mj = floes.mass[i], floes.mass[j]
            c = 2 * bs.damping_ratio * sqrt(k * b.L * mi*mj/(mi + mj))
            Fd = bond_damping_force(floes.α[i], b.attach_self, Ui, floes.ξ[i],
                                    floes.α[j], b.attach_partner, Uj, floes.ξ[j], c)
            τd = cross2(rot(floes.α[i], b.attach_self), Fd)
            bonds.fx  += F[1] + Fd[1]
            bonds.fy  += F[2] + Fd[2]
            bonds.trq += τ + τd
        end
        # coast bonds: spring to a fixed anchor; the coast is rigid, so stiffness matches
        # floe-topography collisions and damping uses the floe's full mass
        # with coastline bonds off, only floe-floe bonds act
        (bs.coastline_bonds && !isempty(bonds.coast)) || continue
        k = consts.E * floes.height[i] / sqrt(floes.area[i]) / bs.stiffness_ratio
        for b in bonds.coast
            F, τ = fixed_bond_force(xi, floes.α[i], b.attach_self, b.anchor, k, b.L)
            c = 2 * bs.damping_ratio * sqrt(k * b.L * floes.mass[i])
            offi = rot(floes.α[i], b.attach_self)
            Fd = -c * (Ui + floes.ξ[i] * perp(offi))
            bonds.fx  += F[1] + Fd[1]
            bonds.fy  += F[2] + Fd[2]
            bonds.trq += τ + cross2(offi, Fd)
        end
    end
    return
end

# domain period in x and y; Inf on axes that aren't periodic (periodic boundaries come in pairs)
function _periodic_lengths(domain, ::Type{FT}) where {FT}
    Lx = domain.east isa PeriodicBoundary ? FT(domain.east.val - domain.west.val) : FT(Inf)
    Ly = domain.north isa PeriodicBoundary ? FT(domain.north.val - domain.south.val) : FT(Inf)
    return Lx, Ly
end

_wrap1(d, L) = isfinite(L) ? d - L*round(d/L) : d      # avoids Inf*0 = NaN
_wrap(d, Lx, Ly) = SA.SVector(_wrap1(d[1], Lx), _wrap1(d[2], Ly))

function prune_dangling_bonds!(floes)
    ids = Set(floes.id)
    for i in eachindex(floes)
        filter!(b -> b.partner_id in ids, floes.bonds[i].floe)
    end
    return
end
