abstract type AbstractOperator end
abstract type AbstractOneElectronOperator <: AbstractOperator end
abstract type AbstractTwoElectronOperator <: AbstractOperator end

struct OverlapOperator <: AbstractOneElectronOperator end
struct KineticOperator <: AbstractOneElectronOperator end
struct NuclearOperator{
		T<:Real,
		V<:AbstractVector{<:Real},
		M<:AbstractMatrix{<:T}
	} <: AbstractOneElectronOperator
	charges::V
	coords::M
end
NuclearOperator(m::AbstractMolecule) = NuclearOperator(m.numbers, m.coords)
Base.length(op::NuclearOperator) = length(op.charges)

struct CoulombOperator <: AbstractTwoElectronOperator end

"""
	create_constant(::AbstractOperator)

Return a multiplication constant for the given operator.

Integral kernels are multiplied by a constant that defaults to one. This is
useful to avoid unnecessary multiplication in a hot loop.
"""
create_constant(::AbstractOperator) = 1
create_constant(::NuclearOperator) = 2π
create_constant(::CoulombOperator) = 2π^2 * sqrt(π)

"""
	create_scratch(operator::AbstractOperator)

Create a scratch tuple suited for the given operator.

Operators require distinct scratches, normally tree to seven
three-dimentional vectors.
"""
create_scratch(operator::AbstractOperator) = (Vector{Float64}(undef, 3),
											  Vector{Float64}(undef, 3),
											  Vector{Float64}(undef, 3))
create_scratch(operator::CoulombOperator) = (Vector{Float64}(undef, 3),
											 Vector{Float64}(undef, 3),
											 Vector{Float64}(undef, 3),
											 Vector{Float64}(undef, 3),
											 Vector{Float64}(undef, 3),
											 Vector{Float64}(undef, 3),
											 Vector{Float64}(undef, 3))
create_scratch(operator::NuclearOperator) = (Vector{Float64}(undef, 3),
											 Vector{Float64}(undef, 3),
											 Vector{Float64}(undef, 3),
											 Vector{Float64}(undef, 3))
