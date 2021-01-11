@time @testset "dist" begin
    @test GaussianBasisFunctions.dist([0, 0, 0], [1, 1, 0]) == sqrt(2)

    carbon_monoxide = Molecule(
        [6, 8],
        [0.000 0.000 0.000;
         0.000 0.000 1.128]',
    )
    basis = build_sto3g(carbon_monoxide)
    @test GaussianBasisFunctions.dist(basis[1], basis[6]) == 1.128
end
