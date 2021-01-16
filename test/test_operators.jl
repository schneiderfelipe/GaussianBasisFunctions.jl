@time @testset "operators" begin
    Ŝ = OverlapOperator()
    @test Ŝ == Ŝ
    @test Ŝ == OverlapOperator()
    @test Ŝ === OverlapOperator()

    K̂ = KineticOperator()
    @test K̂ == K̂
    @test K̂ == KineticOperator()
    @test K̂ === KineticOperator()

    water = Molecule(
        [8, 1, 1],
        [0.00000  0.00000  0.22700;
        0.00000  1.35300 -0.90800;
        0.00000 -1.35300 -0.90800]',
    )
    basis = build_sto3g(water)
    V̂ = NuclearOperator(water)
    @test V̂ == V̂
    @test V̂ == NuclearOperator(water)
    @test V̂ === NuclearOperator(water)
end
