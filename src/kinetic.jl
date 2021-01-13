# TODO: create a generic, numerically integrated, fallback for any callable
# objects
"""
    oei(a::GaussianBasisFunction, b::GaussianBasisFunction, ::KineticOperator)

Compute a kinetic energy matrix element for a pair of basis functions.
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

function _kinetic_kernel(alpha, a_coord, la, ma, na, beta, b_coord, lb, mb, nb)
    l0 = (2 * (lb + mb + nb) + 3) * _overlap_kernel(
        alpha, a_coord, la, ma, na, beta, b_coord, lb, mb, nb
    )
    lplus2 = _overlap_kernel(
        alpha, a_coord, la, ma, na, beta, b_coord, lb + 2, mb, nb
    ) + _overlap_kernel(
        alpha, a_coord, la, ma, na, beta, b_coord, lb, mb + 2, nb
    ) + _overlap_kernel(
        alpha, a_coord, la, ma, na, beta, b_coord, lb, mb, nb + 2
    )
    lminus2 = lb * (lb - one(lb)) * _overlap_kernel(
        alpha, a_coord, la, ma, na, beta, b_coord, lb - 2, mb, nb
    ) + mb * (mb - one(mb)) * _overlap_kernel(
        alpha, a_coord, la, ma, na, beta, b_coord, lb, mb - 2, nb
    ) + nb * (nb - one(nb)) * _overlap_kernel(
        alpha, a_coord, la, ma, na, beta, b_coord, lb, mb, nb - 2
    )

    return beta * (l0 - 2 * beta * lplus2) - lminus2 / 2
end
