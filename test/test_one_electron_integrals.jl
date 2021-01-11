@testset "one_electron_integrals.jl" begin
    α = 1.0

    # Test _compute_primitive_norm_constant
    @test GaussianBasisFunctions._compute_primitive_norm_constant(
        α, 0, 0, 0
    ) == 0.7127054703549902

    # Test _compute_primitive_third_center
    @test GaussianBasisFunctions._compute_primitive_third_center(
        α, [0, 0, 0], 1.0, [0, 0, 1.128]
    ) == [0, 0, 1.128 / 2]

    # Test overlap
    water = Molecule(
        [8, 1, 1],
        [0.00000  0.00000  0.22700;
         0.00000  1.35300 -0.90800;
         0.00000 -1.35300 -0.90800]',
    )
    basis = build_sto3g(water)
    @test overlap(basis) ≈ [1.00000 0.23670 0.00000  0.00000  0.00000  0.05695  0.05695
                            0.23670 1.00000 0.00000  0.00000  0.00000  0.48979  0.48979
                            0.00000 0.00000 1.00000  0.00000  0.00000  0.00000  0.00000
                            0.00000 0.00000 0.00000  1.00000  0.00000  0.30738 -0.30738
                            0.00000 0.00000 0.00000  0.00000  1.00000 -0.25785 -0.25785
                            0.05695 0.48979 0.00000  0.30738 -0.25785  1.00000  0.28279
                            0.05695 0.48979 0.00000 -0.30738 -0.25785  0.28279  1.00000] atol=2e-5
end
