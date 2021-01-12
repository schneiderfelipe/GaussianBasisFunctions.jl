# TODO: create a generic, numerically integrated, fallback for any callable
# objects
"""
    oei(a::GaussianBasisFunction, b::GaussianBasisFunction, ::KineticOperator)

Computes a kinetic energy matrix element for a pair of basis functions.
"""
function oei(a::GaussianBasisFunction, b::GaussianBasisFunction, ::KineticOperator)
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

            term += b.coeffs[j] * N_b * _kinetic_kernel(
                alpha, a.coord, a.l, a.m, a.n, beta, b.coord, b.l, b.m, b.n
            )
        end

        res += a.coeffs[i] * N_a * term
    end

    return res
end

function _kinetic_kernel(alpha, a_coord, a_l, a_m, a_n, beta, b_coord, b_l, b_m, b_n)
    return beta * ((2 * (b_l + b_m + b_n) + 3) * _overlap_kernel(
        alpha, a_coord, a_l, a_m, a_n, beta, b_coord, b_l, b_m, b_n
    ) - 2 * beta * (
        _overlap_kernel(
            alpha, a_coord, a_l, a_m, a_n, beta, b_coord, b_l + 2, b_m, b_n
        ) + _overlap_kernel(
            alpha, a_coord, a_l, a_m, a_n, beta, b_coord, b_l, b_m + 2, b_n
        ) + _overlap_kernel(
            alpha, a_coord, a_l, a_m, a_n, beta, b_coord, b_l, b_m, b_n + 2
        )
    )) - (
        b_l * (b_l - one(b_l)) * _overlap_kernel(
            alpha, a_coord, a_l, a_m, a_n, beta, b_coord, b_l - 2, b_m, b_n
        ) + b_m * (b_m - one(b_m)) * _overlap_kernel(
            alpha, a_coord, a_l, a_m, a_n, beta, b_coord, b_l, b_m - 2, b_n
        ) + b_n * (b_n - one(b_n)) * _overlap_kernel(
            alpha, a_coord, a_l, a_m, a_n, beta, b_coord, b_l, b_m, b_n - 2
        )
    ) / 2
end
