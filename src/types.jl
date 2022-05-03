abstract type AbstractMolecule end

"""
Molecule with atomic numbers and coordinates as columns of a matrix.

In the future this will include charge and multiplicity as well.
"""
# TODO: create a package for input/output of molecular structures
struct Molecule{V, W <: AbstractMatrix{<:AbstractFloat}} <: AbstractMolecule
	numbers::V
	coords::W
end
Base.length(m::Molecule) = length(m.numbers)

# NOTE: We assume molecules are equal if their atomic numbers and coordinates
# are equal. This can lead to surprises when atoms are permuted.
function Base.hash(m::Molecule, h::UInt)
	return hash(m.coords,
		hash(m.numbers,
		hash(:Molecule, h)))
end
function Base.:(==)(m::Molecule, n::Molecule)
	return isequal(m.coords, n.coords) &&
		isequal(m.numbers, n.numbers)
end
abstract type AbstractBasisFunction end

"""
Gaussian basis function with coordinate, alphas, coefficients and angular
momenta.

I might add atomic number and/or orbital name (e.g., "1s") in the future.
It may be benefitial to also separate the angular momenta into integers `l`,
`m` and `n`.
"""
struct GaussianBasisFunction{T<:Real, V<:AbstractVector{<:T}, W<:AbstractVector{<:T}, Z<:AbstractVector{<:T}} <: AbstractBasisFunction
    # NOTE: it is assumed that alphas and coeffs have the same length
    coord::V
    alphas::W
    coeffs::Z  # ::V might become ::W after we substitute build_sto3g with a real basis set parser
    l::Int
    m::Int
    n::Int
end
Base.length(a::GaussianBasisFunction) = length(a.coeffs)

# NOTE: We assume Gaussian basis functions are equal if their data are equal.
# This can lead to surprises if the data in the vectors are permuted.
function Base.hash(a::GaussianBasisFunction, h::UInt)
    return hash(a.n,
        hash(a.m,
        hash(a.l,
        hash(a.coeffs,
        hash(a.alphas,
        hash(a.coord,
        hash(:GaussianBasisFunction, h)))))))
end
function Base.:(==)(a::GaussianBasisFunction, b::GaussianBasisFunction)
    return isequal(a.n, b.n) &&
        isequal(a.m, b.m) &&
        isequal(a.l, b.l) &&
        isequal(a.coeffs, b.coeffs) &&
        isequal(a.alphas, b.alphas) &&
        isequal(a.coord, b.coord)
end
