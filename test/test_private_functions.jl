@time @testset "private functions" begin
    α = 1.0

    # Test _compute_primitive_norm_constant
    @test GaussianBasisFunctions._compute_primitive_norm_constant(
        α, 0, 0, 0
    ) == 0.7127054703549902

    # Test _compute_primitive_third_center
    @test GaussianBasisFunctions._compute_primitive_third_center(
        α, [0, 0, 0], 1.0, [0, 0, 1.128]
    ) == [0, 0, 1.128 / 2]
end
