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
struct GaussianBasisFunction{V, W, I} <: AbstractBasisFunction
	coord::V
	alphas::W
	coeffs::W
	l::I
	m::I
	n::I
end
Base.length(a::GaussianBasisFunction) = length(a.coeffs)
