using GaussianBasisFunctions
using Test

@testset "GaussianBasisFunctions.jl" begin
    include("test_concrete_types.jl")

    include("test_utils.jl")
    include("test_one_electron_integrals.jl")
end
