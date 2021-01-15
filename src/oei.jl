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
    return Hermitian(M)
end

# TODO: create a generic, numerically integrated, fallback for any callable
# objects, not just GaussianBasisFunctions.
"""
    oei(a::GaussianBasisFunction, b::GaussianBasisFunction, operator::AbstractOperator)

Compute a one electron matrix element for a pair of basis functions.
"""
function oei(a::GaussianBasisFunction, b::GaussianBasisFunction, operator::AbstractOperator)
    pa, pb = similar(a.coord), similar(b.coord)

    res = zero(a.coeffs[1])
    for i in 1:length(a)
        alpha = a.alphas[i]
        N_a = _normalization_constant(alpha, a.l, a.m, a.n)

        term = zero(res)
        for j in 1:length(b)
            beta = b.alphas[j]
            N_b = _normalization_constant(beta, b.l, b.m, b.n)

            term += b.coeffs[j] * N_b * integral_kernel!(
                pa, pb, operator, alpha, a.coord, a.l, a.m, a.n, beta, b.coord, b.l, b.m, b.n
            )
        end

        res += a.coeffs[i] * N_a * term
    end

    return res
end
