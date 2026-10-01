# contains rotational , cross2, bond_force, plus loop that accumulates bond force and torques onto floes
@inline rot(α, r::SA.SVector{2}) =
    SA.SVector(cos(α)*r[1] - sin(α)*r[2], sin(α)*r[1] + cos(α)*r[2])

@inline cross2(r::SA.SVector{2}, F::SA.SVector{2}) = r[1]*F[2] - r[2]*F[1]

# force on floe i (floe j gets the negative), plus torques
function bond_force(xi, αi, ri, xj, αj, rj, k, L)
    offi = rot(αi, ri); offj = rot(αj, rj)
    F = -k * L * ((xi + offi) - (xj + offj))
    return F, cross2(offi, F), cross2(offj, -F)
end

# This is for coastline/wall bonds only and acts like a spring attached to a single fixed point. 
function fixed_bond_force(xi, αi, ri, anchor, k, L)
    offi = rot(αi, ri)
    F = -k * L * ((xi + offi) - anchor)
    return F, cross2(offi, F)
end

@inline perp(r::SA.SVector{2}) = SA.SVector(-r[2], r[1])

# damping force on floe i (floe j gets the negative)
function bond_damping_force(αi, ri, Ui, ωi, αj, rj, Uj, ωj, c)
    vi = Ui + ωi * perp(rot(αi, ri))
    vj = Uj + ωj * perp(rot(αj, rj))
    return -c * (vi - vj)
end