# contains rotational , cross2, bond_force, plus loop that accumulates bond force and torques onto floes
@inline rot(α, r::SA.SVector{2}) =
    SA.SVector(cos(α)*r[1] - sin(α)*r[2], sin(α)*r[1] + cos(α)*r[2])

@inline cross2(r::SA.SVector{2}, F::SA.SVector{2}) = r[1]*F[2] - r[2]*F[1]

# Scale the part of `d` along unit normal `n` by `kn` (motion perpendicular to the bond face:
# tension/compression) and the part along the face by `ks` (shear)
@inline function split_scale(d, n, kn, ks)
    dn = dot(d, n) * n
    return kn * dn + ks * (d - dn)
end

#= World-frame unit normal of a floe-floe bond, pointing from floe i toward floe j. It averages
floe i's face normal with the reverse of floe j's, so both floes split the bond into the same
normal and shear parts (and the force stays equal and opposite) even after they rotate
relative to each other. Falls back on floe i's normal if the two normals cancel. =#
@inline function bond_normal(αi, ni, αj, nj)
    m = rot(αi, ni) - rot(αj, nj)
    l = norm(m)
    return l > eps(typeof(l)) ? m / l : rot(αi, ni)
end

# force on floe i (floe j gets the negative), plus torques; ni/nj are the floes' body-frame
# face normals and kn/ks the normal/shear stiffness per unit face length
function bond_force(xi, αi, ri, ni, xj, αj, rj, nj, kn, ks, L)
    offi = rot(αi, ri); offj = rot(αj, rj)
    d = (xi + offi) - (xj + offj)
    F = -L * split_scale(d, bond_normal(αi, ni, αj, nj), kn, ks)
    # With kn ≠ ks, F isn't parallel to d, so it acts at the midpoint of the two anchors to
    # conserve angular momentum
    return F, cross2(offi - d / 2, F), cross2(offj + d / 2, -F)
end

# This is for coastline/wall bonds only and acts like a spring attached to a single fixed point. 
function fixed_bond_force(xi, αi, ri, ni, anchor, kn, ks, L)
    offi = rot(αi, ri)
    F = -L * split_scale((xi + offi) - anchor, rot(αi, ni), kn, ks)
    return F, cross2(offi, F)
end

@inline perp(r::SA.SVector{2}) = SA.SVector(-r[2], r[1])

# damping force on floe i (floe j gets the negative); n is the bond's world-frame unit normal
# and cn/cs the damping coefficients of normal/shear motion
function bond_damping_force(αi, ri, Ui, ωi, αj, rj, Uj, ωj, n, cn, cs)
    vi = Ui + ωi * perp(rot(αi, ri))
    vj = Uj + ωj * perp(rot(αj, rj))
    return -split_scale(vi - vj, n, cn, cs)
end

# Split the force a bond exerts on its floe into a normal stress along the bond's face normal
# (positive in tension, i.e. pulling the floe toward the bond's partner/the coast) and a shear
# stress along the face. Shear is signed so that it has the same value seen from either floe
# of a floe-floe bond. Stresses are force over cross-sectional area: face length L × thickness h.
function bond_stress(F, α, normal_self, L, h)
    n = rot(α, normal_self)
    A = L * h
    return dot(F, n) / A, cross2(n, F) / A
end
