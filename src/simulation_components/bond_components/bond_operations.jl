# add_bond!, remove_bond!, accumulation loop
# must be included AFTER floe, as bond operations rely on the Floe struct

_centroid(f::Floe{FT}) where {FT} = SA.SVector{2,FT}(f.centroid[1], f.centroid[2])

function add_bond!(fi::Floe{FT}, fj::Floe{FT}, p, n, L;
    shift = zero(SA.SVector{2,FT}),
) where {FT}
    # p: world-frame point on the shared face; n: unit normal pointing from fi toward fj
    # shift: offset from fj to the periodic image of fj that shares the face with fi
    any(b -> b.partner_id == fj.id, fi.bonds.floe) &&
        throw(ArgumentError("floes $(fi.id) and $(fj.id) are already bonded"))
    ai = rot(-fi.α, p - _centroid(fi));  aj = rot(-fj.α, p - (_centroid(fj) + shift))
    ni = rot(-fi.α, n);                  nj = rot(-fj.α, -n)
    push!(fi.bonds.floe, FloeBond{FT}(fj.id, ai, aj, ni, nj, L))
    push!(fj.bonds.floe, FloeBond{FT}(fi.id, aj, ai, nj, ni, L))
    return nothing
end

function add_coast_bond!(fi::Floe{FT}, p, n, L) where {FT}
    # p: world-frame point on the face touching the coast (becomes the fixed anchor);
    # n: unit normal pointing from fi toward the coast
    any(b -> b.anchor ≈ p, fi.bonds.coast) &&
        throw(ArgumentError("floe $(fi.id) is already bonded to the coast at $p"))
    ai = rot(-fi.α, p - _centroid(fi))
    ni = rot(-fi.α, n)
    push!(fi.bonds.coast, CoastBond{FT}(ai, p, ni, L))
    return nothing
end

remove_coast_bonds!(fi::Floe) = (empty!(fi.bonds.coast); nothing)

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
        (mb.attach_self == b.attach_partner && mb.attach_partner == b.attach_self &&
            mb.normal_self == b.normal_partner && mb.normal_partner == b.normal_self &&
            mb.L == b.L) || return false
    end
    return true
end

# exterior ring of a floe polygon as closed, counterclockwise list of points
function _ccw_ring(poly, ::Type{FT}) where {FT}
    pts = [SA.SVector{2,FT}(GI.x(p), GI.y(p)) for p in GI.getpoint(GI.getexterior(poly))]
    pts[1] == pts[end] || push!(pts, pts[1])
    sum(cross2(pts[k], pts[k+1]) for k in 1:length(pts)-1) < 0 && reverse!(pts)
    return pts
end

#= Longest face shared by counterclockwise rings ri and rj, i.e. the longest overlap of an
edge of ri with a collinear edge of rj. Returns (midpoint, unit outward normal of ri, length)
or nothing if no overlap is longer than tol (floes touching at a corner aren't bonded). =#
function _shared_face(ri, rj, tol)
    best, bestL = nothing, tol
    for k in 1:length(ri)-1
        a, b = ri[k], ri[k+1]
        ℓ = norm(b - a)
        ℓ ≤ tol && continue
        t = (b - a) / ℓ
        for m in 1:length(rj)-1
            c, e = rj[m], rj[m+1]
            (abs(cross2(t, c - a)) ≤ tol && abs(cross2(t, e - a)) ≤ tol) || continue
            sc, se = dot(c - a, t), dot(e - a, t)
            s0, s1 = max(zero(ℓ), min(sc, se)), min(ℓ, max(sc, se))
            if s1 - s0 > bestL
                bestL = s1 - s0
                # ring is counterclockwise, so the outward normal is the tangent turned clockwise
                best = (a + t * (s0 + s1) / 2, SA.SVector(t[2], -t[1]), bestL)
            end
        end
    end
    return best
end

"""
    bond_adjacent_floes!(floes, domain; tol)

Bond every pair of floes that share a face, e.g. neighboring floes of a tessellated floe
field. Each bond is anchored at the midpoint of the shared face and its length is the face
length. Faces shared across periodic boundaries are bonded to the partner's periodic image.
Faces must overlap by more than `tol` (m) to be bonded, so floes touching only at a corner
are not bonded. Returns the number of bonds added.
"""
function bond_adjacent_floes!(floes::StructArray{<:Floe{FT}}, domain;
    tol = FT(1e-6) * (domain.east.val - domain.west.val),
) where {FT}
    Lx, Ly = _periodic_lengths(domain, FT)
    xshifts = isfinite(Lx) ? (-Lx, zero(FT), Lx) : (zero(FT),)
    yshifts = isfinite(Ly) ? (-Ly, zero(FT), Ly) : (zero(FT),)
    shifts = [SA.SVector{2,FT}(sx, sy) for sx in xshifts for sy in yshifts]
    rings = [_ccw_ring(p, FT) for p in floes.poly]
    lo = [reduce((u, v) -> min.(u, v), r) for r in rings]
    hi = [reduce((u, v) -> max.(u, v), r) for r in rings]
    nadded = 0
    for i in eachindex(floes)
        floes.ghost_id[i] == 0 || continue
        for j in (i + 1):lastindex(floes)
            floes.ghost_id[j] == 0 || continue
            any(b -> b.partner_id == floes.id[j], floes.bonds[i].floe) && continue
            for s in shifts
                # skip images of floe j whose bounding box doesn't touch floe i's
                (all(lo[j] + s .≤ hi[i] .+ tol) && all(lo[i] .≤ hi[j] + s .+ tol)) || continue
                face = _shared_face(rings[i], [q + s for q in rings[j]], tol)
                isnothing(face) && continue
                p, n, L = face
                add_bond!(floes[i], floes[j], p, n, L; shift = s)
                nadded += 1
                break
            end
        end
    end
    return nadded
end

# distance s ≥ -tol along the ray m + s * n to its first crossing of an edge of ring r,
# skipping edge `skip`; Inf if the ray doesn't cross r
function _ray_hit(m, n, r, tol; skip = 0)
    hit = oftype(tol, Inf)
    for k in 1:length(r)-1
        k == skip && continue
        c, w = r[k], r[k+1] - r[k]
        denom = cross2(n, w)
        abs(denom) ≤ eps(norm(w)) && continue      # ray parallel to edge
        s = cross2(c - m, w) / denom
        u = cross2(c - m, n) / denom
        (-tol ≤ s < hit && 0 ≤ u ≤ 1) && (hit = s)
    end
    return hit
end

"""
    bond_floes_to_coast!(floes, domain; max_gap, tol)

Bond floes to the coast (the domain's topography) along every face that touches it. A face
touches the coast if the coast is within `max_gap` (m) of the face's midpoint, looking out
from the face along its outward normal, with no other floe in between, so faces shared
with other floes and faces with another floe between them and the coast aren't bonded.
`max_gap` allows for gaps between the floe field and the coastline, e.g. when they were
outlined from different data. Each face gets its own coast bond, anchored at the face's
midpoint, with the face's length as its length. Faces already bonded are skipped. Faces
must be longer than `tol` (m). Returns the number of bonds added.
"""
function bond_floes_to_coast!(floes::StructArray{<:Floe{FT}}, domain;
    max_gap = zero(FT),
    tol = FT(1e-6) * (domain.east.val - domain.west.val),
) where {FT}
    isempty(domain.topography) && return 0
    coast = [_ccw_ring(p, FT) for p in domain.topography.poly]
    rings = [_ccw_ring(p, FT) for p in floes.poly]
    lo = [reduce((u, v) -> min.(u, v), r) for r in rings]
    hi = [reduce((u, v) -> max.(u, v), r) for r in rings]
    reach = max_gap + tol
    nadded = 0
    for i in eachindex(floes)
        floes.ghost_id[i] == 0 || continue
        ri = rings[i]
        for k in 1:length(ri)-1
            a, b = ri[k], ri[k+1]
            ℓ = norm(b - a)
            ℓ ≤ tol && continue
            t = (b - a) / ℓ
            n = SA.SVector(t[2], -t[1])     # outward normal of a counterclockwise ring
            m = (a + b) / 2
            any(cb -> cb.anchor ≈ m, floes.bonds[i].coast) && continue
            dcoast = minimum(c -> _ray_hit(m, n, c, tol), coast)
            dcoast ≤ reach || continue
            # bounding box of the ray segment from the face to the coast
            e = m + n * dcoast
            rlo, rhi = min.(m, e) .- tol, max.(m, e) .+ tol
            blocked = any(eachindex(floes)) do j
                floes.ghost_id[j] == 0 || return false
                (all(lo[j] .≤ rhi) && all(rlo .≤ hi[j])) || return false
                _ray_hit(m, n, rings[j], tol; skip = j == i ? k : 0) ≤ dcoast
            end
            blocked && continue
            add_coast_bond!(floes[i], m, n, ℓ)
            nadded += 1
        end
    end
    return nadded
end
