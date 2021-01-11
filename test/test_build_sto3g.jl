@time @testset "build_sto3g" begin
    hydrogen = Molecule(
        [1],
        [0 0 0]',
    )
    basisset = build_sto3g(hydrogen)
    @test basisset[1].coord == [0, 0, 0]
    @test basisset[1].alphas == [3.425250914, 0.6239137298, 0.168855404]
    @test basisset[1].coeffs == [0.1543289673, 0.5353281423, 0.4446345422]
    @test basisset[1].l == 0
    @test basisset[1].m == 0
    @test basisset[1].n == 0

    carbon_monoxide = Molecule(
        [6, 8],
        [0.000 0.000 0.000;
         0.000 0.000 1.128]',
    )
    basisset = build_sto3g(carbon_monoxide)
    @test length(basisset) == 10
    @test basisset[1].coord == [0, 0, 0]
    @test basisset[6].coord == [0, 0, 1.128]
end
