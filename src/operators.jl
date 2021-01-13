abstract type AbstractOperator end

struct OverlapOperator <: AbstractOperator end
struct KineticOperator <: AbstractOperator end

struct NuclearOperator{V, W <: AbstractMatrix{<:AbstractFloat}} <: AbstractOperator
	charges::V
	coords::W
end
NuclearOperator(m::AbstractMolecule) = NuclearOperator(m.numbers, m.coords)
Base.length(op::NuclearOperator) = length(op.charges)
