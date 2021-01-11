abstract type AbstractMolecule end

"""
Molecule with atomic numbers and coordinates as columns of a matrix.

In the future this will include charge and multiplicity as well.
"""
# TODO: create a package for input/output of molecular structures
struct Molecule{V, W} <: AbstractMolecule
	numbers::V
	coords::W
end
Base.length(m::Molecule) = length(m.numbers)

abstract type AbstractBasisFunction end

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
