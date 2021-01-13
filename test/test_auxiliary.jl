@time @testset "auxiliary" begin
    α = 1.0

    # Test Boys function
    @test GaussianBasisFunctions.boys.(
        0, [9e-5, 0.1, 1, 4]
    ) ≈ [0.9999700008099822, 0.9676433126355923, 0.7468241328124272, 0.4410406953812108]
    @test GaussianBasisFunctions.boys.(
        1, [0, 9e-5, 1, 4]
    ) ≈ [1/3, 0.3333153339118912, 0.18947234582049235, 0.05284063206155958]

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
