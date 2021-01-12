@time @testset "overlap" begin
    water = Molecule(
        [8, 1, 1],
        [0.00000  0.00000  0.22700;
         0.00000  1.35300 -0.90800;
         0.00000 -1.35300 -0.90800]',
    )
    basis = build_sto3g(water)

    # Tests for pairs of basis functions
    @test overlap(basis[1], basis[end]) ≈ 0.05695 atol=5e-6
    @test overlap(basis[1], basis[1]) == 1.0

    # Test for a complete basis set
    @test overlap(basis) ≈ [1.00000 0.23670 0.00000  0.00000  0.00000  0.05695  0.05695
                            0.23670 1.00000 0.00000  0.00000  0.00000  0.48979  0.48979
                            0.00000 0.00000 1.00000  0.00000  0.00000  0.00000  0.00000
                            0.00000 0.00000 0.00000  1.00000  0.00000  0.30738 -0.30738
                            0.00000 0.00000 0.00000  0.00000  1.00000 -0.25785 -0.25785
                            0.05695 0.48979 0.00000  0.30738 -0.25785  1.00000  0.28279
                            0.05695 0.48979 0.00000 -0.30738 -0.25785  0.28279  1.00000] atol=2e-5

    # Test overlap between primitives
    @test GaussianBasisFunctions._primitive_overlap(
        0.1, [1.0, 2.0, 0.0], 1, 1, 0, 0.2, [2.0, 1.0, 0.0], 0, 1, 0
    ) == 28.55923564396704
end
