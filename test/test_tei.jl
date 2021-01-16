@time @testset "tei" begin
    # Reference tests for a complete basis set
    hydrogen = Molecule(
        [1, 1],
        [0.0 0.0 -0.8;
         0.0 0.0  0.8]',
    )
    basis = build_sto3g(hydrogen)
    G = tei(basis, CoulombOperator())
    @test length(G) == length(basis)^4
    @test G[:, :, 1, 1] ≈ [0.77461 0.38568;
                           0.38568 0.53068] atol=8e-6
    @test G[:, :, 1, 2] ≈ [0.38568 0.23234;
                           0.23234 0.38568] atol=1e-5
    @test G[:, :, 2, 1] ≈ [0.38568 0.23234;
                           0.23234 0.38568] atol=1e-5
    @test G[:, :, 2, 2] ≈ [0.53068 0.38568;
                           0.38568 0.77461] atol=8e-6

    # Reference tests for groups of basis functions
    @test tei(basis[1], basis[2], basis[1], basis[2], CoulombOperator()) ≈ 0.23234 atol=5e-6
end
