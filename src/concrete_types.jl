"""
Molecule with atomic numbers and coordinates as columns of a matrix.
"""
struct Molecule{V, M} <: AbstractMolecule
	numbers::V
	coords::M
end

"""
Gaussian basis function with coordinate, alphas, coefficients and angular
momenta.

I might add atomic number and/or orbital name (e.g., "1s") in the future.
It may be benefitial to also separate the angular momenta into integers `l`,
`m` and `n`.
"""
struct GaussianBasisFunction{A, B} <: AbstractBasisFunction
	coord::A
	alphas::B
	coeffs::B
	anguls::NTuple{3, <:Integer}
end
