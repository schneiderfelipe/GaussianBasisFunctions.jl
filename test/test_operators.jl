@time @testset "operators" begin
    kinetic_op = KineticOperator()
    @test kinetic_op == kinetic_op
    @test kinetic_op == KineticOperator()
    @test kinetic_op === KineticOperator()

    overlap_op = OverlapOperator()
    @test overlap_op == overlap_op
    @test overlap_op == OverlapOperator()
    @test overlap_op === OverlapOperator()
end
