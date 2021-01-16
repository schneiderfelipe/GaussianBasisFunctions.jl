@time @testset "auxiliary" begin
    α = 1.0

    # Test gaussian_norm
    @test GaussianBasisFunctions.gaussian_norm(α, 0, 0, 0) == 0.7127054703549902

    # Test gpt!
    p_coord = Vector{Float64}(undef, 3)
    @test GaussianBasisFunctions.gpt!(p_coord, α, [0, 0, 0], 1.0, [0, 0, 1.128]) == (0.5293041815490781, α + 1)
    @test p_coord == [0, 0, 1.128 / 2]

    # Test c
    @test GaussianBasisFunctions.c(1, 1, 2, 0.1, 0.2) ≈ 0.08

    # Test Sᵢ
    @test GaussianBasisFunctions.Sᵢ(1, 2, 0.1, 0.2, 0.3) == 2.7096468290777316

    # Test _trrt
    @test GaussianBasisFunctions._trrt.([1, 3, 6, 10, 15]) == [1, 2, 3, 4, 5]
end
