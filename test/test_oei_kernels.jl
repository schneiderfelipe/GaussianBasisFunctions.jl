@time @testset "oei kernels" begin
    water = Molecule(
        [8, 1, 1],
        [0.00000  0.00000  0.22700;
        0.00000  1.35300 -0.90800;
        0.00000 -1.35300 -0.90800]',
    )
    basis = build_sto3g(water)

    # Consistency checks for kernels
    @test GaussianBasisFunctions.integral_kernel(
        OverlapOperator(), 0.1, [1.0, 2.0, 0.0], 1, 1, 0, 0.2, [2.0, 1.0, 0.0], 0, 1, 0
    ) == 28.55923564396704
    @test GaussianBasisFunctions.integral_kernel(
        KineticOperator(), 0.1, [1.0, 2.0, 0.0], 1, 1, 0, 0.2, [2.0, 1.0, 0.0], 0, 1, 0
    ) == 12.234093080988277
    @test GaussianBasisFunctions.integral_kernel(
        NuclearOperator(water), 0.1, [1.0, 2.0, 0.0], 1, 1, 0, 0.2, [2.0, 1.0, 0.0], 0, 1, 0
    ) ≈ -48.42642469823712
end
