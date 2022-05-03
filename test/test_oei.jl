@time @testset "oei" begin
    water = Molecule(
        [8, 1, 1],
        [0.00000  0.00000  0.22700;
         0.00000  1.35300 -0.90800;
         0.00000 -1.35300 -0.90800]',
    )
    basis = build_sto3g(water)

    # Reference tests for a complete basis set
    @test oei(basis, OverlapOperator()) ≈ [1.000000 0.236704 0.000000  0.000000  0.000000  0.056954  0.056954;
                                           0.236704 1.000000 0.000000  0.000000  0.000000  0.489792  0.489792;
                                           0.000000 0.000000 1.000000  0.000000  0.000000  0.000000  0.000000;
                                           0.000000 0.000000 0.000000  1.000000  0.000000  0.307379 -0.307379;
                                           0.000000 0.000000 0.000000  0.000000  1.000000 -0.257853 -0.257853;
                                           0.056954 0.489792 0.000000  0.307379 -0.257853  1.000000  0.282791;
                                           0.056954 0.489792 0.000000 -0.307379 -0.257853  0.282791  1.000000] atol=2e-6
    @test oei(basis, KineticOperator()) ≈ [29.003200 -0.168011 0.000000  0.000000  0.000000 -0.000807 -0.000807;
                                           -0.168011  0.808128 0.000000  0.000000  0.000000  0.139790  0.139790;
                                            0.000000  0.000000 2.528731  0.000000  0.000000  0.000000  0.000000;
                                            0.000000  0.000000 0.000000  2.528731  0.000000  0.231386 -0.231386;
                                            0.000000  0.000000 0.000000  0.000000  2.528731 -0.194105 -0.194105;
                                           -0.000807  0.139790 0.000000  0.231386 -0.194105  0.760032  0.016145;
                                           -0.000807  0.139790 0.000000 -0.231386 -0.194105  0.016145  0.760032] atol=5e-6
    @test oei(basis, NuclearOperator(water)) ≈ [-61.750945  -7.451144   0.000000   0.000000   0.020933 -1.846697 -1.846697;
                                                 -7.451144 -10.167129   0.000000   0.000000   0.244043 -4.034052 -4.034052;
                                                  0.000000   0.000000 -10.006114   0.000000   0.000000  0.000000  0.000000;
                                                  0.000000   0.000000   0.000000 -10.160795   0.000000 -2.251996  2.251996;
                                                  0.020933   0.244043   0.000000   0.000000 -10.114965  1.974772  1.974772;
                                                 -1.846697  -4.034052   0.000000  -2.251996   1.974772 -5.944108 -1.832457;
                                                 -1.846697  -4.034052   0.000000   2.251996   1.974772 -1.832457 -5.944108] atol=4e-6

    # Reference tests for pairs of basis functions
    @test oei(basis[1], basis[1], OverlapOperator()) ≈ 1.0
    @test oei(basis[1], basis[end], OverlapOperator()) ≈ 0.056954 atol=5e-7
    @test oei(basis[2], basis[end], KineticOperator()) ≈ 0.139790 atol=5e-7
    @test oei(basis[4], basis[end], NuclearOperator(water)) ≈ 2.251996 atol=2e-7
end
