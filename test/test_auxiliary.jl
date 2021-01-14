@time @testset "auxiliary" begin
    α = 1.0

    # Test _normalization_constant
    @test GaussianBasisFunctions._normalization_constant(
        α, 0, 0, 0,
    ) == 0.7127054703549902

    # Test _third_center
    @test GaussianBasisFunctions._third_center(
        α, [0, 0, 0], 1.0, [0, 0, 1.128], α + 1.0
    ) == [0, 0, 1.128 / 2]

    # Test _ck
    @test GaussianBasisFunctions._ck(
        1, 1, 2, 0.1, 0.2,
    ) ≈ 0.08

    # Test _Si
    @test GaussianBasisFunctions._Si(
        1, 2, 0.1, 0.2, 0.3,
    ) == 2.7096468290777316
end
