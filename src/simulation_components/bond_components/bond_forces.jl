# contains rotational , cross2, bond_force, plus loop that accumulates bond force and torques onto floes
@inline rot(α, r::SVector{2}) =
    SVector(cos(α)*r[1] - sin(α)*r[2], sin(α)*r[1] + cos(α)*r[2])

@inline cross2(r::SVector{2}, F::SVector{2}) = r[1]*F[2] - r[2]*F[1]

# force on floe i (floe j gets the negative), plus torques
function bond_force(xi, αi, ri, xj, αj, rj, k, L)
    offi = rot(αi, ri); offj = rot(αj, rj)
    F = -k * L * ((xi + offi) - (xj + offj))
    return F, cross2(offi, F), cross2(offj, -F)
end