@time @testset "oei kernels" begin
    a_coord = [1.0, 2.0, 0.0]
    b_coord = [2.0, 1.0, 0.0]
    pa, pb = similar(a_coord), similar(b_coord)  # scratch

    # Consistency check for the OverlapOperator kernel
    @test GaussianBasisFunctions.integral_kernel!(
        pa, pb, OverlapOperator(), 0.1, a_coord, 1, 1, 0, 0.2, b_coord, 0, 1, 0
    ) == 28.55923564396704

    # Consistency check for the KineticOperator kernel
    @test GaussianBasisFunctions.integral_kernel!(
        pa, pb, KineticOperator(), 0.1, a_coord, 1, 1, 0, 0.2, b_coord, 0, 1, 0
    ) == 12.234093080988277

    # Consistency check for the NuclearOperator kernel
    water = Molecule(
        [8, 1, 1],
        [0.00000  0.00000  0.22700;
        0.00000  1.35300 -0.90800;
        0.00000 -1.35300 -0.90800]',
    )
    @test GaussianBasisFunctions.integral_kernel!(
        pa, pb, NuclearOperator(water), 0.1, a_coord, 1, 1, 0, 0.2, b_coord, 0, 1, 0
    ) ≈ -48.42642469823712
end
