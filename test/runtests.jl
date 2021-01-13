using GaussianBasisFunctions
using Test

@time @testset "GaussianBasisFunctions.jl" begin
    # Types
    include("test_types.jl")
    include("test_operators.jl")

    # Utilities
    include("test_doublefactorial.jl")
    include("test_dist.jl")
    include("test_build_sto3g.jl")
    include("test_auxiliary.jl")

    # Integrals
    include("test_oei.jl")
    include("test_oei_kernels.jl")
end
