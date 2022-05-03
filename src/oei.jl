# NOTE: This assumes real basis functions, should we have
# `AbstractRealBasisFunction <: AbstractBasisFunction`?
"""
    oei(
        basis::AbstractVector{<:AbstractBasisFunction},
        operator::AbstractOneElectronOperator
    )

Compute a symmetric one electron matrix for a basis set.
"""
function oei(
        basis::AbstractVector{<:AbstractBasisFunction},
        operator::AbstractOneElectronOperator
    )
    n = length(basis)

    # TODO: define eltype for AbstractBasisFunctions so that we avoid typeof
    T = typeof(basis[1].coeffs[1])
    M = Matrix{T}(undef, n, n)

    scratch = create_scratch(operator)
    constant = create_constant(operator)
    @views for j in 1:n, i in 1:j
        M[i, j] = constant * oei!(scratch, basis[i], basis[j], operator)
    end

    # We defined above only the upper triangular part, so we return a wrapped
    # symmetric matrix (uplo == :U is the default).
    return Hermitian(M)
end

"""
    oei(
        a::GaussianBasisFunction,
        b::GaussianBasisFunction,
        operator::AbstractOneElectronOperator
    )

Compute a one electron matrix element for a pair of basis functions.

Please prefer the potentially more efficient `oei!` over this function.
"""
function oei(
        a::GaussianBasisFunction,
        b::GaussianBasisFunction,
        operator::AbstractOneElectronOperator
    )
    return create_constant(operator) * oei!(
        create_scratch(operator), a, b, operator
    )
end

"""
    oei!(
        scratch,
        a::GaussianBasisFunction,
        b::GaussianBasisFunction,
        operator::AbstractOneElectronOperator,
    )

Compute a one electron matrix element for a pair of basis functions.

This function works in-place and overrides `scratch`.
"""
function oei!(
        scratch,
        a::GaussianBasisFunction,
        b::GaussianBasisFunction,
        operator::AbstractOneElectronOperator,
    )
    return combine_gaussians(a, b) do i, j
        return integral_kernel!(scratch, operator, a, b, i, j)
    end
end
