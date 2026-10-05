# add_bond!, remove_bond!, accumulation loop
# must be included AFTER floe, as bond operations rely on the Floe struct

_centroid(f::Floe{FT}) where {FT} = SA.SVector{2,FT}(f.centroid[1], f.centroid[2])

function add_bond!(fi::Floe{FT}, fj::Floe{FT}, p::SA.SVector{2,FT}, L::FT) where {FT}
    any(b -> b.partner_id == fj.id, fi.bonds.floe) &&
        throw(ArgumentError("floes $(fi.id) and $(fj.id) are already bonded"))
    ri = rot(-fi.α, p - _centroid(fi))
    rj = rot(-fj.α, p - _centroid(fj))
    push!(fi.bonds.floe, FloeBond{FT}(fj.id, ri, rj, L))
    push!(fj.bonds.floe, FloeBond{FT}(fi.id, rj, ri, L))
    return nothing
end

function remove_bond!(fi::Floe, fj::Floe)
    filter!(b -> b.partner_id != fj.id, fi.bonds.floe)
    filter!(b -> b.partner_id != fi.id, fj.bonds.floe)
    return nothing
end

# true if every bond on every (non-ghost) floe has an exact mirror on its partner
function bonds_symmetric(floes)
    byid = Dict(f.id => f for f in floes if f.ghost_id == 0)
    for f in values(byid), b in f.bonds.floe
        p = get(byid, b.partner_id, nothing)
        p === nothing && return false
        m = findall(x -> x.partner_id == f.id, p.bonds.floe)
        length(m) == 1 || return false
        mb = p.bonds.floe[only(m)]
        (mb.r_self == b.r_partner && mb.r_partner == b.r_self && mb.L == b.L) || return false
    end
    return true
end