@time @testset "types" begin
    # Test Molecule type
    carbon_monoxide = Molecule(
        [6, 8],
        [0.000 0.000 0.000;
        0.000 0.000 1.128]',
    )
    @test length(carbon_monoxide) == 2
    @test carbon_monoxide == carbon_monoxide
    @test carbon_monoxide == Molecule(
        [6, 8],
        [0.000 0.000 0.000;
        0.000 0.000 1.128]',
    )

    # Test GaussianBasisFunction type
    basisfun = GaussianBasisFunction(
        [0.0, 0.0, 0.0],
        [3.425250914, 0.6239137298, 0.168855404],
        [0.1543289673, 0.5353281423, 0.4446345422],
        0, 0, 0,
    )
    @test length(basisfun) == 3
    @test basisfun == basisfun
    @test basisfun == GaussianBasisFunction(
        [0.0, 0.0, 0.0],
        [3.425250914, 0.6239137298, 0.168855404],
        [0.1543289673, 0.5353281423, 0.4446345422],
        0, 0, 0,
    )
end
