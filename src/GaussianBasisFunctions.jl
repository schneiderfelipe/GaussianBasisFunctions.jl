module GaussianBasisFunctions

using LinearAlgebra: Hermitian
using LinearAlgebra: hermitian
using LinearAlgebra: hermitian_type
using LinearAlgebra: checksquare
using LinearAlgebra: char_uplo
using LinearAlgebra: sym_uplo
using SpecialFunctions: erf
using SpecialFunctions: gamma
using SpecialFunctions: gamma_inc

export Molecule
export GaussianBasisFunction
export GTensor
export OverlapOperator
export KineticOperator
export NuclearOperator

export build_sto3g

export boys
export oei

# Types
include("types.jl")
include("operators.jl")

# Utilities
include("utils.jl")
include("build_sto3g.jl")
include("auxiliary.jl")
include("boys.jl")
include("gtensor.jl")

# Integrals
include("oei.jl")
include("oei_kernels.jl")

end
