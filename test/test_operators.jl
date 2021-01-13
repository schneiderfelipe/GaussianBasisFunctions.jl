@time @testset "operators" begin
    overlap_op = OverlapOperator()
    @test overlap_op == overlap_op
    @test overlap_op == OverlapOperator()
    @test overlap_op === OverlapOperator()

    kinetic_op = KineticOperator()
    @test kinetic_op == kinetic_op
    @test kinetic_op == KineticOperator()
    @test kinetic_op === KineticOperator()

    water = Molecule(
        [8, 1, 1],
        [0.00000  0.00000  0.22700;
        0.00000  1.35300 -0.90800;
        0.00000 -1.35300 -0.90800]',
    )
    basis = build_sto3g(water)
    nuclear_op = NuclearOperator(water)
    @test nuclear_op == nuclear_op
    @test nuclear_op == NuclearOperator(water)
    @test nuclear_op === NuclearOperator(water)
end
