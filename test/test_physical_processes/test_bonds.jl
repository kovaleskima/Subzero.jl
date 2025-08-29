@testset "Bonds" begin
    @testset "Fracture Criteria" begin
        FT = Float64
        # Test BondFractures criteria
        @test BondFractures() isa BondFractures

    end
end