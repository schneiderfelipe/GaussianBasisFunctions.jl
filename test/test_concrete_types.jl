@testset "concrete_types.jl" begin
    # Test Base.length
    basis = GaussianBasisFunction(
        [0, 0, 0],
        [3.425250914, 0.6239137298, 0.168855404],
        [0.1543289673, 0.5353281423, 0.4446345422],
        0, 0, 0,
    )
    @test length(basis) == 3
end
