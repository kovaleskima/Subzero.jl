# test/test_simulation_components/bond_components/test_bonds.jl
using Test, Subzero, StaticArrays, StructArrays
using LinearAlgebra: norm, dot
using Subzero: add_bond!, remove_bond!, bonds_symmetric, rot, add_coast_bond!, remove_coast_bonds!

@testset "bond bookkeeping" begin
    A = Floe{Float64}([[[0.0,0.0],[0.0,5.0],[10.0,5.0],[10.0,0.0],[0.0,0.0]]], 1.0; id = 1)
    B = Floe{Float64}([[[10.0,0.0],[10.0,5.0],[20.0,5.0],[20.0,0.0],[10.0,0.0]]], 1.0; id = 2)
    p = SVector(10.0, 2.5)
    n = SVector(1.0, 0.0)   # unit normal from A toward B

    @testset "add_bond! creates mirrored bonds in body frames" begin
        add_bond!(A, B, p, n, 5.0)
        @test A.bonds.floe[1].partner_id == 2 && B.bonds.floe[1].partner_id == 1
        @test A.bonds.floe[1].attach_self ≈ SVector(5.0, 0.0)
        @test B.bonds.floe[1].attach_self ≈ SVector(-5.0, 0.0)
        @test norm(A.bonds.floe[1].normal_self) ≈ 1
        @test bonds_symmetric([A, B])
    end

    @testset "duplicate bond is rejected" begin
        @test_throws ArgumentError add_bond!(A, B, p, n, 5.0)
    end

    @testset "bond normals" begin
        # reuses the A–B bond created above
        nA = A.bonds.floe[1].normal_self
        nB = B.bonds.floe[1].normal_self

        @test norm(nA) ≈ 1
        @test norm(nB) ≈ 1

        # world normals are opposite right after creation
        @test rot(A.α, nA) ≈ -rot(B.α, nB)

        # normal is perpendicular to the face (here the face is vertical, tangent = (0, 1))
        @test abs(dot(nA, SVector(0.0, 1.0))) < 1e-12

        # points out of A, toward B (convex floes: away from A's centroid)
        @test dot(rot(A.α, nA), p - SVector(A.centroid[1], A.centroid[2])) > 0
    end

    @testset "body-frame vector round-trips through rotation" begin
        A2 = Floe{Float64}([[[0.0,0.0],[0.0,5.0],[10.0,5.0],[10.0,0.0],[0.0,0.0]]], 1.0; id = 3, α = π/2)
        B2 = Floe{Float64}([[[10.0,0.0],[10.0,5.0],[20.0,5.0],[20.0,0.0],[10.0,0.0]]], 1.0; id = 4)
        add_bond!(A2, B2, p, n, 5.0)
        world = SVector(A2.centroid[1], A2.centroid[2]) + rot(A2.α, A2.bonds.floe[1].attach_self)
        @test world ≈ p
    end

    @testset "remove_bond! clears both sides" begin
        remove_bond!(A, B)
        @test isempty(A.bonds) && isempty(B.bonds)
        @test bonds_symmetric([A, B])
    end

    @testset "symmetry check catches a broken mirror" begin
        add_bond!(A, B, p, n, 5.0)
        empty!(B.bonds.floe)                 # simulate a one-sided removal bug
        @test !bonds_symmetric([A, B])
    end
end

@testset "coast bond bookkeeping" begin
    # floe C sits right of a wall along x = 0
    C = Floe{Float64}([[[0.0,0.0],[0.0,5.0],[10.0,5.0],[10.0,0.0],[0.0,0.0]]], 1.0; id = 5)
    p = SVector(0.0, 2.5)
    n = SVector(-1.0, 0.0)   # unit normal from C toward the wall

    @testset "add_coast_bond! stores a body-frame attachment and a fixed world anchor" begin
        add_coast_bond!(C, p, n, 5.0)
        b = only(C.bonds.coast)
        @test b.attach_self ≈ SVector(-5.0, 0.0)
        @test b.anchor == p
        @test b.normal_self ≈ n
        @test norm(b.normal_self) ≈ 1
        @test b.L == 5.0
        @test isempty(C.bonds.floe) && !isempty(C.bonds)
    end

    @testset "duplicate coast bond is rejected" begin
        @test_throws ArgumentError add_coast_bond!(C, p, n, 5.0)
    end

    @testset "remove_coast_bonds! clears coast bonds" begin
        remove_coast_bonds!(C)
        @test isempty(C.bonds)
    end
end

@testset "coastline_bonds turns coast bond forces on and off" begin
    grid = RegRectilinearGrid(; x0 = -1e2, xf = 1e2, y0 = -1e2, yf = 1e2, Δx = 10, Δy = 10)
    domain = Subzero.Domain(;
        north = OpenBoundary(North; grid), south = OpenBoundary(South; grid),
        east = OpenBoundary(East; grid), west = OpenBoundary(West; grid),
    )
    consts = Constants()
    # A–B bonded to each other, C bonded to a wall along x = 30
    A = Floe{Float64}([[[0.0,0.0],[0.0,5.0],[10.0,5.0],[10.0,0.0],[0.0,0.0]]], 1.0; id = 1)
    B = Floe{Float64}([[[10.0,0.0],[10.0,5.0],[20.0,5.0],[20.0,0.0],[10.0,0.0]]], 1.0; id = 2)
    C = Floe{Float64}([[[30.0,0.0],[30.0,5.0],[40.0,5.0],[40.0,0.0],[30.0,0.0]]], 1.0; id = 3)
    add_bond!(A, B, SVector(10.0, 2.5), SVector(1.0, 0.0), 5.0)
    add_coast_bond!(C, SVector(30.0, 2.5), SVector(-1.0, 0.0), 5.0)
    B.centroid[1] += 1.0     # stretch the floe-floe bond
    C.centroid[1] += 1.0     # stretch the coast bond
    floes = StructArray([A, B, C])

    Subzero.timestep_bonds!(floes, domain, consts, BondSettings(bonds_on = true), 10)
    @test floes.bonds[1].fx > 0 && floes.bonds[2].fx < 0    # floe-floe bond still acts
    @test floes.bonds[3].fx == 0 && floes.bonds[3].fy == 0 && floes.bonds[3].trq == 0

    Subzero.timestep_bonds!(floes, domain, consts, BondSettings(bonds_on = true, coastline_bonds = true), 10)
    @test floes.bonds[1].fx > 0 && floes.bonds[2].fx < 0
    @test floes.bonds[3].fx < 0                              # pulled back toward the wall
end

@testset "bond_adjacent_floes!" begin
    using Subzero: bond_adjacent_floes!
    rect(x0, y0, w, h) = [[[x0, y0], [x0, y0 + h], [x0 + w, y0 + h], [x0 + w, y0], [x0, y0]]]
    floes_from(coords) = StructArray([Floe{Float64}(c, 1.0; id = k) for (k, c) in enumerate(coords)])
    ext = (; x0 = 0.0, xf = 100.0, y0 = 0.0, yf = 100.0)
    collision_domain = Subzero.Domain(;
        north = CollisionBoundary(North; ext...), south = CollisionBoundary(South; ext...),
        east = CollisionBoundary(East; ext...), west = CollisionBoundary(West; ext...),
    )
    periodic_domain = Subzero.Domain(;
        north = CollisionBoundary(North; ext...), south = CollisionBoundary(South; ext...),
        east = PeriodicBoundary(East; ext...), west = PeriodicBoundary(West; ext...),
    )

    @testset "shared face is bonded at its midpoint" begin
        floes = floes_from([rect(10, 10, 10, 5), rect(20, 10, 10, 5)])
        @test bond_adjacent_floes!(floes, collision_domain) == 1
        b = floes.bonds[1].floe[1]
        @test b.partner_id == 2 && b.L ≈ 5
        @test b.attach_self ≈ SVector(5.0, 0.0) && b.attach_partner ≈ SVector(-5.0, 0.0)
        @test b.normal_self ≈ SVector(1.0, 0.0)
        @test bonds_symmetric(floes)
    end

    @testset "bond length is the overlap of partially shared faces" begin
        floes = floes_from([rect(10, 10, 10, 10), rect(20, 15, 10, 10)])
        @test bond_adjacent_floes!(floes, collision_domain) == 1
        b = floes.bonds[1].floe[1]
        @test b.L ≈ 5
        @test SVector(floes.centroid[1]...) + b.attach_self ≈ SVector(20.0, 17.5)
    end

    @testset "floes touching at a corner or apart aren't bonded" begin
        floes = floes_from([rect(10, 10, 10, 10), rect(20, 20, 10, 10), rect(50, 50, 10, 10)])
        @test bond_adjacent_floes!(floes, collision_domain) == 0
        @test all(isempty, floes.bonds)
    end

    @testset "grid of floes and existing bonds" begin
        floes = floes_from([rect(10i, 10j, 10, 10) for i in 1:3 for j in 1:3])
        @test bond_adjacent_floes!(floes, collision_domain) == 12
        @test bonds_symmetric(floes)
        # running again doesn't duplicate bonds
        @test bond_adjacent_floes!(floes, collision_domain) == 0
    end

    @testset "faces shared across periodic boundaries" begin
        floes = floes_from([rect(0, 10, 10, 10), rect(90, 10, 10, 10)])
        @test bond_adjacent_floes!(floes, collision_domain) == 0
        @test bond_adjacent_floes!(floes, periodic_domain) == 1
        b = floes.bonds[1].floe[1]
        @test b.L ≈ 10 && b.normal_self ≈ SVector(-1.0, 0.0)
        # partner anchor sits on the face shared with floe 1, at the partner's periodic image
        @test SVector(floes.centroid[1]...) + b.attach_self ≈ SVector(0.0, 15.0)
        @test SVector(floes.centroid[2]...) + SVector(-100.0, 0.0) + b.attach_partner ≈ SVector(0.0, 15.0)
        @test bonds_symmetric(floes)
    end
end
