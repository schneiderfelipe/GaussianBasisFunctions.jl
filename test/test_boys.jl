@time @testset "boys" begin
    Δx = 2e-3
    x = 1e-8:Δx:20
    for n in 0:4
        @test maximum(abs.(boys.(n, x) - GaussianBasisFunctions._boys_gamma.(n, x))) ≤ √eps()
    end
end
