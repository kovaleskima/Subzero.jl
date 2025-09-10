@testset "Bonds" begin
    @testset "Fracture Criteria" begin
        FT = Float64
        # Test BondFractures criteria
        @test BondFractures() isa BondFractures
    end
    @testset "Initialization" begin
        FT = Float64

        @test NoBond() isa NoBond
    end
end