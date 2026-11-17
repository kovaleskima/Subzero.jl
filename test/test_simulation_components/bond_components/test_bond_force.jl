# test/test_bond_force.jl
using Test
using StaticArrays
using Random
using LinearAlgebra: dot, norm
using Subzero: bond_force, bond_normal, fixed_bond_force, rot, cross2

const k, L = 2.5e3, 40.0      # arbitrary but nonzero, so scaling errors show up
const kL = k * L
const kn, ks = 4.0e3, 1.0e3   # distinct normal and shear stiffness

@testset "bond_force" begin

    # Reference config: bond center coincides on both floes (unstretched)
    xi, αi, ri = SVector(0.0, 0.0), 0.3, SVector(1.0, 0.5)
    xj, αj = SVector(3.0, 1.0), -0.2
    rj = rot(-αj, (xi + rot(αi, ri)) - xj)        # chosen so p_i == p_j exactly
    # face normals of both floes in their body frames, for a world normal n pointing i -> j
    n = SVector(cos(0.4), sin(0.4))
    ni, nj = rot(-αi, n), rot(-αj, -n)
    t = SVector(-n[2], n[1])                      # world tangent along the face
    # with equal normal and shear stiffness, the bond is an isotropic spring
    iso(xi, αi, ri, xj, αj, rj) = bond_force(xi, αi, ri, ni, xj, αj, rj, nj, k, k, L)

    @testset "unstretched bond gives zero force and torque" begin
        for (a, b) in ((k, k), (kn, ks))
            F, τi, τj = bond_force(xi, αi, ri, ni, xj, αj, rj, nj, a, b, L)
            @test F ≈ SVector(0.0, 0.0) atol = 1e-12 * kL
            @test τi ≈ 0 atol = 1e-12 * kL
            @test τj ≈ 0 atol = 1e-12 * kL
        end
    end

    @testset "translating floe j pulls floe i toward it, equal and opposite" begin
        δ = SVector(0.1, -0.2)
        F, τi, τj = iso(xi, αi, ri, xj + δ, αj, rj)
        @test F ≈ k * L * δ
        # swapping the roles must give the negative force
        Fswap, _, _ = bond_force(xj + δ, αj, rj, nj, xi, αi, ri, ni, k, k, L)
        @test Fswap ≈ -F
    end

    @testset "normal and shear stiffness act on their own parts of the stretch" begin
        a, b = 0.1, -0.3
        F, _, _ = bond_force(xi, αi, ri, ni, xj + a * n, αj, rj, nj, kn, ks, L)
        @test F ≈ kn * L * a * n                  # pulled apart across the face
        F, _, _ = bond_force(xi, αi, ri, ni, xj + b * t, αj, rj, nj, kn, ks, L)
        @test F ≈ ks * L * b * t                  # slid along the face
        F, _, _ = bond_force(xi, αi, ri, ni, xj + a * n + b * t, αj, rj, nj, kn, ks, L)
        @test dot(F, n) ≈ kn * L * a && dot(F, t) ≈ ks * L * b
        # equal and opposite from either floe
        Fswap, _, _ = bond_force(xj + a * n + b * t, αj, rj, nj, xi, αi, ri, ni, kn, ks, L)
        @test Fswap ≈ -F
    end

    @testset "bond normal is shared by both floes after relative rotation" begin
        @test bond_normal(αi, ni, αj, nj) ≈ n
        Δ = 0.2   # rotate floe j relative to floe i: normal sits halfway between the two
        nij = bond_normal(αi, ni, αj + Δ, nj)
        @test norm(nij) ≈ 1
        @test nij ≈ rot(Δ / 2, n)
        @test bond_normal(αj + Δ, nj, αi, ni) ≈ -nij
    end

    @testset "rigid rotation of both floes about a common point" begin
        c = SVector(5.0, -2.0)
        θ = 0.7
        Rθ(v) = rot(θ, v)
        xi′ = c + Rθ(xi - c);  αi′ = αi + θ
        xj′ = c + Rθ(xj - c);  αj′ = αj + θ

        # unstretched stays unstretched
        F, _, _ = iso(xi′, αi′, ri, xj′, αj′, rj)
        @test F ≈ SVector(0.0, 0.0) atol = 1e-10 * kL

        # stretched: the force rotates covariantly with the frame, isotropic or not
        δ = SVector(0.3, 0.1)
        for (a, b) in ((k, k), (kn, ks))
            F0, _, _ = bond_force(xi, αi, ri, ni, xj + δ, αj, rj, nj, a, b, L)
            F1, _, _ = bond_force(xi′, αi′, ri, ni, c + Rθ(xj + δ - c), αj′, rj, nj, a, b, L)
            @test F1 ≈ Rθ(F0)
        end
    end

    @testset "rotating floe i alone about its centroid" begin
        Δα = 1e-3
        ri0 = SVector(1.0, 0.0)
        xi0, xj0, rj0 = SVector(0.0, 0.0), SVector(2.0, 0.0), SVector(-1.0, 0.0)
        n0 = SVector(1.0, 0.0)
        F, τi, τj = bond_force(xi0, Δα, ri0, n0, xj0, 0.0, rj0, -n0, k, k, L)
        # exact: p_i - p_j = (cos Δα - 1, sin Δα)
        @test F ≈ -k * L * SVector(cos(Δα) - 1, sin(Δα))
        # first order: stretch is perpendicular to r_i with size Δα * |r_i|
        @test F ≈ -k * L * Δα * SVector(0.0, 1.0) rtol = 1e-2
        # that stretch is along the face, so it is resisted by the shear stiffness
        F, _, _ = bond_force(xi0, Δα, ri0, n0, xj0, 0.0, rj0, -n0, kn, ks, L)
        @test F ≈ -ks * L * Δα * SVector(0.0, 1.0) rtol = 1e-2
    end

    @testset "total torque about the origin vanishes (random configs)" begin
        rng = MersenneTwister(1)
        for _ in 1:100
            local xi, xj, ri, rj, αi, αj, ni, nj    # don't overwrite the reference config
            xi, xj = SVector(randn(rng, 2)...), SVector(randn(rng, 2)...)
            ri, rj = SVector(randn(rng, 2)...), SVector(randn(rng, 2)...)
            αi, αj = 2π * rand(rng), 2π * rand(rng)
            θ = 2π * rand(rng)
            ni, nj = SVector(cos(θ), sin(θ)), SVector(-cos(θ), -sin(θ))
            F, τi, τj = bond_force(xi, αi, ri, ni, xj, αj, rj, nj, kn, ks, L)
            total = τi + τj + cross2(xi, F) + cross2(xj, -F)
            @test abs(total) < 1e-8 * max(1, kn * L)
        end
    end

    @testset "fixed (coast) bond splits along the floe's face normal" begin
        a, b = 0.1, -0.3
        anchor = xi + rot(αi, ri)
        F, _ = fixed_bond_force(xi, αi, ri, ni, anchor - a * n - b * t, kn, ks, L)
        @test dot(F, n) ≈ -kn * L * a && dot(F, t) ≈ -ks * L * b
    end

    @testset "type stability and Float32 support" begin
        args = (SVector(0f0, 0f0), 0.1f0, SVector(1f0, 0f0), SVector(1f0, 0f0),
                SVector(2f0, 0f0), 0.2f0, SVector(-1f0, 0f0), SVector(-1f0, 0f0), 2f3, 5f2, 40f0)
        @inferred bond_force(args...)
        F, τi, τj = bond_force(args...)
        @test eltype(F) === Float32
        @test τi isa Float32
    end
end
