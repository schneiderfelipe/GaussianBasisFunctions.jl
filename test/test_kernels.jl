@time @testset "oei kernels" begin
    a = GaussianBasisFunction(
        [1.0, 2.0, 0.0],
        [0.1],
        [1.0],
        1, 1, 0,
    )
    b = GaussianBasisFunction(
        [2.0, 1.0, 0.0],
        [0.2],
        [1.0],
        0, 1, 0,
    )

    # Consistency check for the OverlapOperator kernel
    Ŝ = OverlapOperator()
    @test GaussianBasisFunctions.integral_kernel!(GaussianBasisFunctions.create_scratch(Ŝ), Ŝ, a, b, 1, 1) == 28.55923564396704

    # Consistency check for the KineticOperator kernel
    K̂ = KineticOperator()
    @test GaussianBasisFunctions.integral_kernel!(GaussianBasisFunctions.create_scratch(K̂), K̂, a, b, 1, 1) == 12.234093080988277

    # Consistency check for the CoulombOperator kernel
    Ĝ = CoulombOperator()
    @test GaussianBasisFunctions.integral_kernel!(GaussianBasisFunctions.create_scratch(Ĝ), Ĝ, a, a, b, b, 1, 1, 1, 1) ≈ 0

    # Consistency check for the NuclearOperator kernel
    water = Molecule(
        [8, 1, 1],
        [0.00000  0.00000  0.22700;
        0.00000  1.35300 -0.90800;
        0.00000 -1.35300 -0.90800]',
    )
    V̂ = NuclearOperator(water)
    @test GaussianBasisFunctions.integral_kernel!(GaussianBasisFunctions.create_scratch(V̂), V̂, a, b, 1, 1) ≈ -48.42642469823712 / 2π
end
