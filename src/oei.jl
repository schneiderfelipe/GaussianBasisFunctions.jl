# NOTE: This assumes real basis functions, should we have
# `AbstractRealBasisFunction <: AbstractBasisFunction`?
"""
    oei(basis::AbstractVector{<:AbstractBasisFunction}, operator::AbstractOperator)

Compute a symmetric one electron matrix for a basis set.
"""
function oei(basis::AbstractVector{<:AbstractBasisFunction}, operator::AbstractOperator)
    n = length(basis)

    # TODO: define eltype for AbstractBasisFunctions so that we avoid typeof
    T = typeof(basis[1].coeffs[1])
    M = Matrix{T}(undef, n, n)

    for j in 1:n, i in 1:j
        M[i, j] = oei(basis[i], basis[j], operator)
    end

    # We defined above only the lower triangular part, so we return a wrapped
    # symmetric matrix.
    return Symmetric(M)
end
