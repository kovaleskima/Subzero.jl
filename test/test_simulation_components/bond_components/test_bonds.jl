# test/test_simulation_components/bond_components/test_bonds.jl
using Test, Subzero, StaticArrays
using Subzero: add_bond!, remove_bond!, bonds_symmetric, rot

@testset "bond bookkeeping" begin
    A = Floe{Float64}([[[0.0,0.0],[0.0,5.0],[10.0,5.0],[10.0,0.0],[0.0,0.0]]], 1.0; id = 1)
    B = Floe{Float64}([[[10.0,0.0],[10.0,5.0],[20.0,5.0],[20.0,0.0],[10.0,0.0]]], 1.0; id = 2)
    p = SVector(10.0, 2.5)

    @testset "add_bond! creates mirrored bonds in body frames" begin
        add_bond!(A, B, p, 5.0)
        @test A.bonds.floe[1].partner_id == 2 && B.bonds.floe[1].partner_id == 1
        @test A.bonds.floe[1].r_self ≈ SVector(5.0, 0.0)
        @test B.bonds.floe[1].r_self ≈ SVector(-5.0, 0.0)
        @test bonds_symmetric([A, B])
    end

    @testset "duplicate bond is rejected" begin
        @test_throws ArgumentError add_bond!(A, B, p, 5.0)
    end

    @testset "body-frame vector round-trips through rotation" begin
        A2 = Floe{Float64}([[[0.0,0.0],[0.0,5.0],[10.0,5.0],[10.0,0.0],[0.0,0.0]]], 1.0; id = 3, α = π/2)
        B2 = Floe{Float64}([[[10.0,0.0],[10.0,5.0],[20.0,5.0],[20.0,0.0],[10.0,0.0]]], 1.0; id = 4)
        add_bond!(A2, B2, p, 5.0)
        world = SVector(A2.centroid[1], A2.centroid[2]) + rot(A2.α, A2.bonds.floe[1].r_self)
        @test world ≈ p
    end

    @testset "remove_bond! clears both sides" begin
        remove_bond!(A, B)
        @test isempty(A.bonds) && isempty(B.bonds)
        @test bonds_symmetric([A, B])
    end

    @testset "symmetry check catches a broken mirror" begin
        add_bond!(A, B, p, 5.0)
        empty!(B.bonds.floe)                 # simulate a one-sided removal bug
        @test !bonds_symmetric([A, B])
    end
end