# test/test_bonds.jl
using Test, Subzero
using StaticArrays
using Random

const k, L = 2.5e3, 40.0      # arbitrary but nonzero, so scaling errors show up

@testset "bond_force" begin

    # Reference config: bond center coincides on both floes (unstretched)
    xi, αi, ri = SVector(0.0, 0.0), 0.3, SVector(1.0, 0.5)
    xj, αj = SVector(3.0, 1.0), -0.2
    rj = rot(-αj, (xi + rot(αi, ri)) - xj)        # chosen so p_i == p_j exactly

    @testset "unstretched bond gives zero force and torque" begin
        F, τi, τj = bond_force(xi, αi, ri, xj, αj, rj, k, L)
        @test F ≈ SVector(0.0, 0.0) atol = 1e-12
        @test τi ≈ 0 atol = 1e-12
        @test τj ≈ 0 atol = 1e-12
    end

    @testset "translating floe j pulls floe i toward it, equal and opposite" begin
        δ = SVector(0.1, -0.2)
        F, τi, τj = bond_force(xi, αi, ri, xj + δ, αj, rj, k, L)
        @test F ≈ k * L * δ
        # swapping the roles must give the negative force
        Fswap, _, _ = bond_force(xj + δ, αj, rj, xi, αi, ri, k, L)
        @test Fswap ≈ -F
    end

    @testset "rigid rotation of both floes about a common point" begin
        c = SVector(5.0, -2.0)
        θ = 0.7
        Rθ(v) = rot(θ, v)
        xi′ = c + Rθ(xi - c);  αi′ = αi + θ
        xj′ = c + Rθ(xj - c);  αj′ = αj + θ

        # unstretched stays unstretched
        F, _, _ = bond_force(xi′, αi′, ri, xj′, αj′, rj, k, L)
        @test F ≈ SVector(0.0, 0.0) atol = 1e-10

        # stretched: the force rotates covariantly with the frame
        δ = SVector(0.3, 0.1)
        F0, _, _ = bond_force(xi, αi, ri, xj + δ, αj, rj, k, L)
        F1, _, _ = bond_force(xi′, αi′, ri, c + Rθ(xj + δ - c), αj′, rj, k, L)
        @test F1 ≈ Rθ(F0)
    end

    @testset "rotating floe i alone about its centroid" begin
        Δα = 1e-3
        ri0 = SVector(1.0, 0.0)
        xi0, xj0, rj0 = SVector(0.0, 0.0), SVector(2.0, 0.0), SVector(-1.0, 0.0)
        F, τi, τj = bond_force(xi0, Δα, ri0, xj0, 0.0, rj0, k, L)
        # exact: p_i - p_j = (cos Δα - 1, sin Δα)
        @test F ≈ -k * L * SVector(cos(Δα) - 1, sin(Δα))
        # first order: stretch is perpendicular to r_i with size Δα * |r_i|
        @test F ≈ -k * L * Δα * SVector(0.0, 1.0) rtol = 1e-2
    end

    @testset "total torque about the origin vanishes (random configs)" begin
        rng = MersenneTwister(1)
        for _ in 1:100
            xi, xj = SVector(randn(rng, 2)...), SVector(randn(rng, 2)...)
            ri, rj = SVector(randn(rng, 2)...), SVector(randn(rng, 2)...)
            αi, αj = 2π * rand(rng), 2π * rand(rng)
            F, τi, τj = bond_force(xi, αi, ri, xj, αj, rj, k, L)
            total = τi + τj + cross2(xi, F) + cross2(xj, -F)
            @test abs(total) < 1e-8 * max(1, k * L)
        end
    end

    @testset "type stability and Float32 support" begin
        args = (SVector(0f0, 0f0), 0.1f0, SVector(1f0, 0f0),
                SVector(2f0, 0f0), 0.2f0, SVector(-1f0, 0f0), 2f3, 40f0)
        @inferred bond_force(args...)
        F, τi, τj = bond_force(args...)
        @test eltype(F) === Float32
        @test τi isa Float32
    end
end