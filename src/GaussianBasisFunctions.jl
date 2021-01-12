module GaussianBasisFunctions

using LinearAlgebra: Symmetric

export Molecule
export GaussianBasisFunction

export build_sto3g
export overlap
export kinetic

include("types.jl")

include("utils.jl")

include("build_sto3g.jl")

include("auxiliary.jl")
include("overlap.jl")
include("kinetic.jl")

end
