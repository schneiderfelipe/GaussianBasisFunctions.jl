using GaussianBasisFunctions
using Test

@time @testset "GaussianBasisFunctions.jl" begin
    include("test_types.jl")

    include("test_doublefactorial.jl")
    include("test_dist.jl")

    include("test_build_sto3g.jl")

    include("test_auxiliary.jl")
    include("test_overlap.jl")
end
