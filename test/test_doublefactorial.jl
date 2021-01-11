@time @testset "doublefactorial" begin
    @test GaussianBasisFunctions.doublefactorial(5) == 15
    @test GaussianBasisFunctions.doublefactorial(6) == 48
    @test GaussianBasisFunctions.doublefactorial(7) == 105
end
