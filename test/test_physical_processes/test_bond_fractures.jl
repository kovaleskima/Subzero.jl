using Test, Subzero, StaticArrays, StructArrays
using Subzero: add_bond!, add_coast_bond!, bonds_symmetric

@testset "Bond Fractures" begin
    @testset "MohrCoulomb envelope" begin
        env = MohrCoulomb(; cohesion = 1e5, μ = 0.5, tensile_strength = 5e4)
        fails(σn, τ) = Subzero._mohr_coulomb_fails(env, σn, τ)
        @test !fails(0.0, 0.0)
        @test !fails(0.0, 0.99e5) && fails(0.0, 1.01e5)    # cohesion sets shear strength
        @test !fails(-1e5, 1.4e5)                          # compression strengthens
        @test fails(4e4, 0.9e5)                            # tension weakens
        @test !fails(4.9e4, 0.0) && fails(5.1e4, 0.0)      # tension cutoff
        # default tensile strength is the envelope's apex
        @test MohrCoulomb(; cohesion = 1e5, μ = 0.5).tensile_strength == 2e5
        @test MohrCoulomb(; μ = 0).tensile_strength == Inf
        @test MohrCoulomb(Float32).cohesion isa Float32
        @test_throws ArgumentError MohrCoulomb(; cohesion = -1)
        @test_throws ArgumentError MohrCoulomb(; μ = -1)
        @test_throws ArgumentError MohrCoulomb(; tensile_strength = -1)
        # coast envelope defaults to the floe envelope
        crit = MohrCoulombBondCriteria(; floe = env)
        @test crit.coast === env
    end

    grid = RegRectilinearGrid(; x0 = -1e2, xf = 1e2, y0 = -1e2, yf = 1e2, Δx = 10, Δy = 10)
    domain = Subzero.Domain(;
        north = OpenBoundary(North; grid), south = OpenBoundary(South; grid),
        east = OpenBoundary(East; grid), west = OpenBoundary(West; grid),
    )
    consts = Constants()
    bs = BondSettings(; bonds_on = true, coastline_bonds = true)
    # A–B bonded across x = 10, C bonded to a wall along x = 30, all 1 m thick
    function bonded_floes(; δB = SVector(0.0, 0.0), δC = SVector(0.0, 0.0))
        A = Floe{Float64}([[[0.0,0.0],[0.0,5.0],[10.0,5.0],[10.0,0.0],[0.0,0.0]]], 1.0; id = 1)
        B = Floe{Float64}([[[10.0,0.0],[10.0,5.0],[20.0,5.0],[20.0,0.0],[10.0,0.0]]], 1.0; id = 2)
        C = Floe{Float64}([[[30.0,0.0],[30.0,5.0],[40.0,5.0],[40.0,0.0],[30.0,0.0]]], 1.0; id = 3)
        add_bond!(A, B, SVector(10.0, 2.5), SVector(1.0, 0.0), 5.0)
        add_coast_bond!(C, SVector(30.0, 2.5), SVector(-1.0, 0.0), 5.0)
        B.centroid .+= δB
        C.centroid .+= δC
        return StructArray([A, B, C])
    end
    # with no velocities there is no damping, so σ = k L δ / (L h) = k δ
    k_floe = consts.E / (2sqrt(50.0)) / bs.stiffness_ratio
    k_coast = consts.E / sqrt(50.0) / bs.stiffness_ratio
    δ = 0.01

    @testset "bond stress" begin
        floes = bonded_floes(; δB = SVector(δ, 0.0))
        b = floes.bonds[1].floe[1]
        xi = SVector(floes.centroid[1]...)
        F, _ = Subzero._floe_bond_force(floes, 1, 2, b, xi, SVector(0.0, 0.0), Inf, Inf, consts, bs)
        σn, τ = Subzero.bond_stress(F, floes.α[1], b.normal_self, b.L, 1.0)
        @test σn ≈ k_floe * δ                          # pulled apart -> tension
        @test τ ≈ 0 atol = 1e-8
        # same stress seen from the partner floe
        b2 = floes.bonds[2].floe[1]
        F2, _ = Subzero._floe_bond_force(floes, 2, 1, b2, SVector(floes.centroid[2]...), SVector(0.0, 0.0), Inf, Inf, consts, bs)
        @test all(Subzero.bond_stress(F2, floes.α[2], b2.normal_self, b2.L, 1.0) .≈ (σn, τ))
        # sliding along the face -> pure shear
        floes = bonded_floes(; δB = SVector(0.0, δ))
        b = floes.bonds[1].floe[1]
        F, _ = Subzero._floe_bond_force(floes, 1, 2, b, xi, SVector(0.0, 0.0), Inf, Inf, consts, bs)
        σn, τ = Subzero.bond_stress(F, floes.α[1], b.normal_self, b.L, 1.0)
        @test σn ≈ 0 atol = 1e-8
        @test abs(τ) ≈ k_floe * δ
        # coast bond pulled away from the wall -> tension
        floes = bonded_floes(; δC = SVector(δ, 0.0))
        cb = floes.bonds[3].coast[1]
        F, _ = Subzero._coast_bond_force(floes, 3, cb, SVector(floes.centroid[3]...), SVector(0.0, 0.0), consts, bs)
        σn, τ = Subzero.bond_stress(F, floes.α[3], cb.normal_self, cb.L, 1.0)
        @test σn ≈ k_coast * δ
        @test τ ≈ 0 atol = 1e-8
    end

    @testset "stress accumulation" begin
        floes = bonded_floes(; δB = SVector(δ, 0.0))
        b = floes.bonds[1].floe[1]
        acc(b, bs) = Subzero._accumulated_bond_stress(b, bs)
        @test acc(b, bs) == (0.0, 0.0)                       # new bond carries no stress
        # timestep_bonds! records the same stress on both sides of a floe-floe bond
        Subzero.timestep_bonds!(floes, domain, consts, bs, 10)
        @test only(b.stress_history)[1] ≈ k_floe * δ
        @test only(floes.bonds[2].floe[1].stress_history) ≈ only(b.stress_history)
        # one timestep only counts for its share of the full window's weights
        r = exp(-bs.stress_decay)
        wsum = sum(r^k for k in 0:bs.stress_nhistory - 1)
        @test acc(b, bs)[1] ≈ k_floe * δ / wsum
        # history is capped and reaches the full stress under sustained load
        for _ in 1:15
            Subzero.timestep_bonds!(floes, domain, consts, bs, 10)
        end
        @test length(b.stress_history) == bs.stress_nhistory
        @test acc(b, bs)[1] ≈ k_floe * δ
        # recent timesteps get more weight: newest is double -> pulled above the plain mean
        hist = pushfirst!(fill(SVector(1.0, 1.0), 9), SVector(2.0, -1.0))
        b_hist = Subzero.FloeBond{Float64}(2, b.attach_self, b.attach_partner, b.normal_self, b.L, hist)
        weights = [r^k for k in 0:9]
        @test acc(b_hist, bs)[1] ≈ sum(weights .* first.(hist)) / wsum
        @test acc(b_hist, bs)[1] > 1.1
        @test acc(b_hist, bs)[2] ≈ abs(sum(weights .* last.(hist))) / wsum   # signed shear averaged
        bs_flat = BondSettings(; bonds_on = true, stress_decay = 0)
        @test acc(b_hist, bs_flat)[1] ≈ 1.1
        # invalid accumulation settings fall back to defaults
        bad = @test_logs (:warn,) (:warn,) BondSettings(; stress_nhistory = 0, stress_decay = -1)
        @test bad.stress_nhistory == 10 && bad.stress_decay == 0.2
    end

    @testset "fracture_bonds!" begin
        σ_floe, σ_coast = k_floe * δ, k_coast * δ
        weak(σ) = MohrCoulomb(; cohesion = 10σ, μ = 0.5, tensile_strength = 0.9σ)
        strong(σ) = MohrCoulomb(; cohesion = 10σ, μ = 0.5, tensile_strength = 1.1σ)
        # hold floes in place and record `n` timesteps of bond stress
        function loaded(n; δB = SVector(δ, 0.0), δC = SVector(δ, 0.0), bs = bs)
            floes = bonded_floes(; δB, δC)
            for _ in 1:n
                Subzero.timestep_bonds!(floes, domain, consts, bs, 10)
            end
            return floes
        end
        n = bs.stress_nhistory

        # NoBondFracture never breaks bonds
        floes = loaded(n)
        @test Subzero.fracture_bonds!(floes, NoBondFracture(), bs) == 0
        @test length(floes.bonds[1].floe) == 1 && length(floes.bonds[3].coast) == 1

        # one timestep of stress isn't enough to break a bond, a sustained load is
        crit = MohrCoulombBondCriteria(; floe = weak(σ_floe), coast = weak(σ_coast))
        @test Subzero.fracture_bonds!(loaded(1), crit, bs) == 0
        @test Subzero.fracture_bonds!(loaded(n), crit, bs) == 2

        # floe and coast bonds use their own envelopes
        floes = loaded(n)
        crit = MohrCoulombBondCriteria(; floe = weak(σ_floe), coast = strong(σ_coast))
        @test Subzero.fracture_bonds!(floes, crit, bs) == 1
        @test isempty(floes.bonds[1].floe) && isempty(floes.bonds[2].floe)
        @test bonds_symmetric(floes)
        @test length(floes.bonds[3].coast) == 1

        floes = loaded(n)
        crit = MohrCoulombBondCriteria(; floe = strong(σ_floe), coast = weak(σ_coast))
        @test Subzero.fracture_bonds!(floes, crit, bs) == 1
        @test length(floes.bonds[1].floe) == 1 && length(floes.bonds[2].floe) == 1
        @test isempty(floes.bonds[3].coast)

        # coast bonds aren't loaded or checked when coastline bonds are off
        bs_no_coast = BondSettings(; bonds_on = true)
        floes = loaded(n; bs = bs_no_coast)
        @test isempty(floes.bonds[3].coast[1].stress_history)
        @test Subzero.fracture_bonds!(floes, MohrCoulombBondCriteria(; floe = strong(σ_floe), coast = weak(σ_coast)),
            bs_no_coast) == 0
        @test length(floes.bonds[3].coast) == 1

        # compression doesn't break a bond with a tiny tensile strength
        floes = loaded(n; δB = SVector(-δ, 0.0), δC = SVector(0.0, 0.0))
        crit = MohrCoulombBondCriteria(; floe = MohrCoulomb(; cohesion = σ_floe, μ = 0.5, tensile_strength = 0.0))
        @test Subzero.fracture_bonds!(floes, crit, bs) == 0
        # but shear above the cohesion does
        floes = loaded(n; δB = SVector(0.0, δ), δC = SVector(0.0, 0.0))
        crit = MohrCoulombBondCriteria(; floe = MohrCoulomb(; cohesion = 0.9σ_floe, μ = 0.5))
        @test Subzero.fracture_bonds!(floes, crit, bs) == 1
        @test bonds_symmetric(floes)
    end
end
