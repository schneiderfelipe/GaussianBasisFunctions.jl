module GaussianBasisFunctions

using LinearAlgebra: Symmetric

export Molecule
export GaussianBasisFunction

export build_sto3g
export overlap

include("types.jl")

include("utils.jl")

include("auxiliary.jl")
include("overlap.jl")

end
