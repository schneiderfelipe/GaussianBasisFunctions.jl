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
export CoulombOperator
export NuclearOperator

export build_sto3g

export boys
export oei
export tei

# Types
include("types.jl")
include("operators.jl")

# Utilities
include("utils.jl")
include("build_sto3g.jl")
include("auxiliary.jl")
include("boys.jl")
include("gtensor.jl")
include("loops.jl")

# Integrals
include("oei.jl")
include("tei.jl")
include("kernels.jl")

end
