# TODO: create a generic, numerically integrated, fallback for any callable
# objects
# TODO: this might benefic from merging with oei(..., ::KineticOperator) and
# calling a kernel.
"""
    oei(a::GaussianBasisFunction, b::GaussianBasisFunction, ::OverlapOperator)

Compute an overlap matrix element for a pair of basis functions.
"""
function oei(a::GaussianBasisFunction, b::GaussianBasisFunction, ::OverlapOperator)
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

            term += b.coeffs[j] * N_b * _overlap_kernel(
                alpha, a.coord, a.l, a.m, a.n, beta, b.coord, b.l, b.m, b.n
            )
        end

        res += a.coeffs[i] * N_a * term
    end

    return res
end

"""
    _overlap_kernel(alpha, a_coord, la, ma, na, beta, b_coord, lb, mb, nb)

Compute the overlap between two primitive Gaussian basis functions.
"""
@inline function _overlap_kernel(alpha, a_coord, la, ma, na, beta, b_coord, lb, mb, nb)
    # TODO: should we compare inputs and return one if they are the same?

    gamma = alpha + beta

    # TODO: the three lines below are bottlenecks and might benefit from
    # inplace operations (in order to avoid allocations). This could be
    # accomplished by reusing intermediate vectors by the caller.
    p_coord = _third_center(alpha, a_coord, beta, b_coord, gamma)
    pa = @. p_coord - a_coord
    pb = @. p_coord - b_coord

    K_p = _gaussian_prod_factor(alpha, a_coord, beta, b_coord, gamma)

    S_x = _Si(la, lb, pa[1], pb[1], gamma)
    S_y = _Si(ma, mb, pa[2], pb[2], gamma)
    S_z = _Si(na, nb, pa[3], pb[3], gamma)

    return K_p * S_x * S_y * S_z
end
