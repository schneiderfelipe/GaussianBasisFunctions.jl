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
    include("test_boys.jl")
    include("test_gtensor.jl")

    # Integrals
    include("test_oei.jl")
    include("test_tei.jl")
    include("test_kernels.jl")
end
