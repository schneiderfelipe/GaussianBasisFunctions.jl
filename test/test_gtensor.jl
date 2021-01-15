@time @testset "gtensor" begin
    G = GTensor{Float64}(undef, 2)  # uplo == :U by default
    H = GTensor{Float64}(undef, 2, :L)

    @test size(G) == (2, 2, 2, 2)
    @test size(H) == size(G)
    for d in 1:4
        @test size(G, d) == 2
        @test size(H, d) == size(G, d)
    end
    @test size(G, 5) == 1
    @test size(H, 5) == size(G, 5)

    # Set all six unique elements in lexicographic order.
    G[1, 1, 1, 1] = 1
    G[1, 1, 1, 2] = 2
    G[1, 1, 2, 2] = 3
    G[1, 2, 1, 2] = 4
    G[1, 2, 2, 2] = 5
    G[2, 2, 2, 2] = 6

    H[1, 1, 1, 1] = 1
    H[1, 1, 1, 2] = 2
    H[1, 1, 2, 2] = 3
    H[1, 2, 1, 2] = 4
    H[1, 2, 2, 2] = 5
    H[2, 2, 2, 2] = 6

    # Test all three unique slices.
    @test G[:, :, 1, 1] == [1 2; 2 3]
    @test G[:, :, 2, 1] == [2 4; 4 5]
    @test G[:, :, 2, 2] == [3 5; 5 6]

    @test H[:, :, 1, 1] == G[:, :, 1, 1]
    @test H[:, :, 2, 1] == G[:, :, 2, 1]
    @test H[:, :, 2, 2] == G[:, :, 2, 2]

    # Test the only symmetry not yet tested.
    @test G[:, :, 2, 1] == G[:, :, 1, 2]
    @test H[:, :, 2, 1] == H[:, :, 1, 2]

    # Assert that the underlying matrices are as expected.
    @test G.data[1, :] == [1, 2, 3]
    @test G.data[2, 2:end] == [4, 5]
    @test G.data[3, 3] == 6

    @test H.data[:, 1] == G.data[1, :]
    @test H.data[2:end, 2] == G.data[2, 2:end]
    @test H.data[3, 3] == G.data[3, 3]

    # Test _compact_indices for symmetric 2 × 2 indices.
    @test GaussianBasisFunctions._compact_indices('U', 2, 1, 1) == 1
    @test GaussianBasisFunctions._compact_indices('U', 2, 1, 2) == 2
    @test GaussianBasisFunctions._compact_indices('U', 2, 2, 1) == 2
    @test GaussianBasisFunctions._compact_indices('U', 2, 2, 2) == 3

    @test GaussianBasisFunctions._compact_indices('L', 2, 1, 1) == 1
    @test GaussianBasisFunctions._compact_indices('L', 2, 2, 1) == 2
    @test GaussianBasisFunctions._compact_indices('L', 2, 1, 2) == 2
    @test GaussianBasisFunctions._compact_indices('L', 2, 2, 2) == 3

    # Test _compact_indices for symmetric 3 × 3 indices.
    @test GaussianBasisFunctions._compact_indices('U', 3, 1, 1) == 1
    @test GaussianBasisFunctions._compact_indices('U', 3, 1, 2) == 2
    @test GaussianBasisFunctions._compact_indices('U', 3, 2, 1) == 2
    @test GaussianBasisFunctions._compact_indices('U', 3, 2, 2) == 3
    @test GaussianBasisFunctions._compact_indices('U', 3, 1, 3) == 4
    @test GaussianBasisFunctions._compact_indices('U', 3, 3, 1) == 4
    @test GaussianBasisFunctions._compact_indices('U', 3, 2, 3) == 5
    @test GaussianBasisFunctions._compact_indices('U', 3, 3, 2) == 5
    @test GaussianBasisFunctions._compact_indices('U', 3, 3, 3) == 6

    @test GaussianBasisFunctions._compact_indices('L', 3, 1, 1) == 1
    @test GaussianBasisFunctions._compact_indices('L', 3, 2, 1) == 2
    @test GaussianBasisFunctions._compact_indices('L', 3, 3, 1) == 3
    @test GaussianBasisFunctions._compact_indices('L', 3, 1, 2) == 2
    @test GaussianBasisFunctions._compact_indices('L', 3, 2, 2) == 4
    @test GaussianBasisFunctions._compact_indices('L', 3, 3, 2) == 5
    @test GaussianBasisFunctions._compact_indices('L', 3, 1, 3) == 3
    @test GaussianBasisFunctions._compact_indices('L', 3, 2, 3) == 5
    @test GaussianBasisFunctions._compact_indices('L', 3, 3, 3) == 6
end
