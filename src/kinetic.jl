# TODO: this is code repetition since `overlap(basis)` is almost the same.
# Should we define `one_electron_operator` and singletons
# `KineticOperator()`, etc.?
"""
    kinetic(basis::AbstractVector{<:AbstractBasisFunction})

Computes a kinetic energy matrix for a basis set.
"""
function kinetic(basis::AbstractVector{<:AbstractBasisFunction})
    n = length(basis)

    T = typeof(basis[1].coeffs[1])
    K = Matrix{T}(undef, n, n)

    for j in 1:n, i in 1:j
        # NOTE: bottleneck of this function
        K[i, j] = kinetic(basis[i], basis[j])
    end

    # We defined the lower triangular part, so we return a symmetric matrix
    return Symmetric(K)
end

# TODO: create a generic, numerically integrated, fallback for any callable
# objects
"""
    kinetic(a::GaussianBasisFunction, b::GaussianBasisFunction)

Computes a kinetic energy matrix element for a pair of basis functions.
"""
function kinetic(a::GaussianBasisFunction, b::GaussianBasisFunction)
    # TODO: can we speed up if invert a and b based on the number of
    # primitives?

    res = zero(a.coeffs[1])
    for i in 1:length(a)
        alpha = a.alphas[i]
        N_a = _normalization_constant(alpha, a.l, a.m, a.n)

        term = zero(res)
        for j in 1:length(b)
            beta = b.alphas[j]
            N_b = _normalization_constant(beta, b.l, b.m, b.n)

            # NOTE: bottleneck of this function
            kern = beta * ((2 * (b.l + b.m + b.n) + 3) * _primitive_overlap(
                alpha, a.coord, a.l, a.m, a.n, beta, b.coord, b.l, b.m, b.n
            ) - 2 * beta * (
                _primitive_overlap(
                    alpha, a.coord, a.l, a.m, a.n, beta, b.coord, b.l + 2, b.m, b.n
                ) + _primitive_overlap(
                    alpha, a.coord, a.l, a.m, a.n, beta, b.coord, b.l, b.m + 2, b.n
                ) + _primitive_overlap(
                    alpha, a.coord, a.l, a.m, a.n, beta, b.coord, b.l, b.m, b.n + 2
                )
            )) - (
                b.l * (b.l - one(b.l)) * _primitive_overlap(
                    alpha, a.coord, a.l, a.m, a.n, beta, b.coord, b.l - 2, b.m, b.n
                ) + b.m * (b.m - one(b.m)) * _primitive_overlap(
                    alpha, a.coord, a.l, a.m, a.n, beta, b.coord, b.l, b.m - 2, b.n
                ) + b.n * (b.n - one(b.n)) * _primitive_overlap(
                    alpha, a.coord, a.l, a.m, a.n, beta, b.coord, b.l, b.m, b.n - 2
                )
            ) / 2

            term += b.coeffs[j] * N_b * kern
        end

        res += a.coeffs[i] * N_a * term
    end

    return res
end
