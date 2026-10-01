# test/test_simulation_components/bond_components/test_bonds.jl
using Test
using StaticArrays
using Subzero: BondSet

# three bonds in a chain: floes 1-2, 2-3, 3-4
function make_bonds(::Type{FT} = Float64) where {FT}
    BondSet(
        [1, 2, 3],                                   # i
        [2, 3, 4],                                   # j
        [SVector{2,FT}(1, 0), SVector{2,FT}(0, 1), SVector{2,FT}(-1, 0)],  # ri
        [SVector{2,FT}(-1, 0), SVector{2,FT}(0, -1), SVector{2,FT}(1, 0)], # rj
        FT[10, 20, 30],                              # L
        [true, true, true],                          # alive
    )
end

@testset "BondSet" begin

    @testset "construction keeps all fields aligned" begin
        b = make_bonds()
        n = length(b.i)
        @test n == 3
        @test length(b.j) == length(b.ri) == length(b.rj) == length(b.L) == length(b.alive) == n
        @test b.i == [1, 2, 3]
        @test b.j == [2, 3, 4]
        @test b.L == [10.0, 20.0, 30.0]
        @test all(b.alive)
    end

    @testset "float type is carried through" begin
        b32 = make_bonds(Float32)
        @test b32 isa BondSet{Float32}
        @test eltype(b32.L) === Float32
        @test eltype(eltype(b32.ri)) === Float32
        @test make_bonds() isa BondSet{Float64}
    end

    @testset "alive flag marks bonds dead without removing them" begin
        b = make_bonds()
        b.alive[2] = false
        @test length(b.i) == 3                  # nothing deleted
        @test count(b.alive) == 2
        @test findall(b.alive) == [1, 3]
        # indices and IDs of the other bonds are unchanged
        @test (b.i[3], b.j[3]) == (3, 4)
    end

    @testset "bonds reference distinct floes, one bond per pair" begin
        b = make_bonds()
        @test all(b.i .!= b.j)                  # no self-bonds
        pairs = [minmax(b.i[n], b.j[n]) for n in eachindex(b.i)]
        @test length(unique(pairs)) == length(pairs)   # no duplicate pair
    end

    # --- optional: only if you add an inner constructor that validates input ---
    # @testset "mismatched lengths are rejected" begin
    #     @test_throws ArgumentError BondSet([1, 2], [2], [SVector(0.0, 0.0)],
    #                                        [SVector(0.0, 0.0)], [1.0], [true])
    # end
end