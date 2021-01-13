module GaussianBasisFunctions

using LinearAlgebra: Symmetric
using SpecialFunctions: gamma
using SpecialFunctions: gamma_inc

export Molecule
export GaussianBasisFunction
export OverlapOperator
export KineticOperator
export NuclearOperator

export build_sto3g

export oei

# Types
include("types.jl")
include("operators.jl")

# Utilities
include("utils.jl")
include("build_sto3g.jl")
include("auxiliary.jl")

# Integrals
include("oei.jl")
include("oei_kernels.jl")

end
