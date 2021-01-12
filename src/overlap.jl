# TODO: create a generic, numerically integrated, fallback for any callable
# objects
# TODO: this might benefic from merging with oei(..., ::KineticOperator) and
# calling a kernel.
"""
    oei(a::GaussianBasisFunction, b::GaussianBasisFunction, ::OverlapOperator)

Computes an overlap matrix element for a pair of basis functions.
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
    _overlap_kernel(alpha, a_coord, a_l, a_m, a_n, beta, b_coord, b_l, b_m, b_n)

Compute the overlap between two primitive Gaussian basis functions.
"""
@inline function _overlap_kernel(alpha, a_coord, a_l, a_m, a_n, beta, b_coord, b_l, b_m, b_n)
    # TODO: should we compare inputs and return one if they are the same?

    gamma = alpha + beta

    # TODO: the three lines below are bottlenecks and might benefit from
    # inplace operations (in order to avoid allocations). This could be
    # accomplished by reusing intermediate vectors by the caller.
    p_coord = _third_center(alpha, a_coord, beta, b_coord, gamma)
    pa = @. p_coord - a_coord
    pb = @. p_coord - b_coord

    S_x = _Si(a_l, b_l, pa[1], pb[1], gamma)
    S_y = _Si(a_m, b_m, pa[2], pb[2], gamma)
    S_z = _Si(a_n, b_n, pa[3], pb[3], gamma)

    # TODO: the following has an extra unneeded ^2
    K_p = exp(-(alpha * beta / gamma) * dist(a_coord, b_coord)^2)

    return K_p * S_x * S_y * S_z
end
