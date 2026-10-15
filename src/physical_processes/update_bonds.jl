# total force (spring + damping) that floe-floe bond `b` exerts on floe i, and its torque
# about floe i's centroid; floe j is the bond partner, xi/Ui are floe i's position/velocity
function _floe_bond_force(floes, i, j, b, xi, Ui, Lx, Ly, consts, bs)
    FT = eltype(xi)
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
    return F + Fd, τ + τd
end

# total force (spring + damping) that coast bond `b` exerts on floe i, and its torque about
# floe i's centroid; the coast is rigid, so stiffness matches floe-topography collisions and
# damping uses the floe's full mass
function _coast_bond_force(floes, i, b, xi, Ui, consts, bs)
    k = consts.E * floes.height[i] / sqrt(floes.area[i]) / bs.stiffness_ratio
    F, τ = fixed_bond_force(xi, floes.α[i], b.attach_self, b.anchor, k, b.L)
    c = 2 * bs.damping_ratio * sqrt(k * b.L * floes.mass[i])
    offi = rot(floes.α[i], b.attach_self)
    Fd = -c * (Ui + floes.ξ[i] * perp(offi))
    return F + Fd, τ + cross2(offi, Fd)
end

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
            F, τ = _floe_bond_force(floes, i, j, b, xi, Ui, Lx, Ly, consts, bs)
            h = min(floes.height[i], floes.height[j])   # the thinner floe is the weak link
            _record_bond_stress!(b, bond_stress(F, floes.α[i], b.normal_self, b.L, h), bs)
            bonds.fx  += F[1]
            bonds.fy  += F[2]
            bonds.trq += τ
        end
        # coast bonds: spring to a fixed anchor
        # with coastline bonds off, only floe-floe bonds act
        (bs.coastline_bonds && !isempty(bonds.coast)) || continue
        for b in bonds.coast
            F, τ = _coast_bond_force(floes, i, b, xi, Ui, consts, bs)
            _record_bond_stress!(b, bond_stress(F, floes.α[i], b.normal_self, b.L, floes.height[i]), bs)
            bonds.fx  += F[1]
            bonds.fy  += F[2]
            bonds.trq += τ
        end
    end
    return
end

# push this timestep's (σn, τ) onto the front of the bond's history, keeping the newest
# `stress_nhistory` timesteps
function _record_bond_stress!(b, (σn, τ), bs)
    h = b.stress_history
    pushfirst!(h, eltype(h)(σn, τ))
    length(h) > bs.stress_nhistory && resize!(h, bs.stress_nhistory)
    return
end

#= Exponentially weighted average of a bond's recorded stress, where the stress from k
timesteps ago has weight exp(-stress_decay * k). Timesteps before the bond existed count as
zero stress, so a new bond needs sustained stress to reach its full accumulated value. Shear is
averaged signed and then made a magnitude, so shear that flips direction cancels out. =#
function _accumulated_bond_stress(b, bs)
    FT = typeof(bs.stress_decay)
    r = exp(-bs.stress_decay)
    σ = zero(SA.SVector{2,FT})
    w = one(FT)
    for s in b.stress_history
        σ += w * s
        w *= r
    end
    wsum = bs.stress_decay == 0 ? FT(bs.stress_nhistory) : (1 - r^bs.stress_nhistory) / (1 - r)
    return σ[1] / wsum, abs(σ[2]) / wsum
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
